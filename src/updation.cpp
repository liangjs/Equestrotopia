#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <cassert>
#include <list>
#include <vector>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <thread>
#include <chrono>
#include <mutex>

#include "light.h"
#include "input.h"
#include "kdtree.h"
#include "updation.h"

#include "lodepng.h"

using namespace std;
using namespace Equestria;

ptnKDTree *ptntree = NULL;
list<string> filelist;
vector<Hitpoint> hits;
vector<Photon *> ptns;
int totforks = 0;
int fd2[2];
mutex filelistlock;
bool endflag = false;
int Nemit = 0;
char cwd[1000];
Camera camera;

void pushFile(const char *str)
{
    filelistlock.lock();
    //printf("pushFile %s\n", str);
    filelist.push_back(str);
    filelistlock.unlock();
}

void listenTerminal()
{
    char line[MAXNLINE];
    int n;
    FILE *file = fdopen(fd2[0], "r");
    while (fscanf(file, "%s", line) != EOF)
        pushFile(line);
    printf("listenTerminal done\n");
    endflag = true;
}

void loadHitpoints()
{
    printf("loading hitpoints...\n");
    FILE *file = fopen("hitpoints.txt", "rb");
    if (file == NULL) {
        fprintf(stderr, "cannot open hitpoints.txt\n");
        exit(1);
    }
    int n;
    fread(&n, 4, 1, file);
    void *buf = malloc(148 * n), *ptr = buf;
    if (fread(buf, 148, n, file) != n) {
        fprintf(stderr, "hitpoints error!!!\n");
        exit(1);
    }
    fclose(file);
    for (int i = 0; i < n; ++i) {
        Hitpoint cur;
        cur.position.Read(ptr);
        cur.normv.Read(ptr);
        cur.raydir.Read(ptr);
        memcpy(&cur.mtl_label, ptr, 4), ptr = (char *)ptr + 4;
        memcpy(&cur.u, ptr, 8), ptr = (char*)ptr + 8;
        memcpy(&cur.v, ptr, 8), ptr = (char*)ptr + 8;
        memcpy(&cur.x, ptr, 4), ptr = (char *)ptr + 4;
        memcpy(&cur.y, ptr, 4), ptr = (char *)ptr + 4;
        cur.wgt.Read(ptr);
        cur.direct.Read(ptr);
        cur.radius = INITRADIUS;
        hits.push_back(cur);
    }
    free(buf);
    printf("hitpoints loaded\n");
}

void clearptnlist()
{
    while (!ptns.empty()) {
        delete ptns.back();
        ptns.pop_back();
    }
}

void readfile()
{
    filelistlock.lock();
    clearptnlist();
    const char *curfile = filelist.front().c_str();
    FILE *file = fopen(curfile, "r");
    int N;
    fread(&N, 4, 1, file);
    Nemit += N;
    filelist.pop_front();
    Photon *curptn = new Photon;
    while (curptn->light.bgn.Read(file) != 0
           && curptn->light.vec.Read(file) != 0
           && curptn->rgb.Read(file) != 0
           && fread(&curptn->mtl, 4, 1, file) != 0) {
        ptns.push_back(curptn);
        curptn = new Photon;
    }
    delete curptn;
    fclose(file);
    filelistlock.unlock();
}

int cntHitUpdated[UPDATE_THREADS];
typedef vector<Hitpoint>::iterator vHit_it;
void __updateHit(const vHit_it &bg, const vHit_it &ed, int id)
{
    for (auto i = bg; i != ed; ++i) {
        auto &ht = *i;
        vector<Photon *> nearptns;
        ptntree->find(ht, nearptns);
        if (!nearptns.empty()) {
            int M = nearptns.size();
            double N = ht.ptncount;
            ht.ptncount += ALPHA * M;
            double ratio = ht.ptncount / (N + M);
            ht.radius *= sqrt(ratio);
            auto &mt = material[ht.mtl_label];
            Point Kd = mt.mtl.getKd(ht.u, ht.v);
            for (auto &ptn : nearptns)
                if (dotsProduct(ptn->light.vec, ht.normv) > EPS) {
                    Point brdf = mt.BRDF(ptn->light.vec, ht.raydir, ht.normv);
                    if (mt.mtl.mapKd != -1)
                        brdf = (brdf + Kd) / 2;
                    ht.tau += elemMult(ptn->rgb, brdf) * dotsProduct(ptn->light.vec, ht.normv);
                }
            ht.tau *= ratio;
        }
        ++cntHitUpdated[id];
    }
};
void updateHitpoints()
{
    ptntree = new ptnKDTree(ptns.begin(), ptns.end());

    if (hits.size() < UPDATE_THREADS)
        fprintf(stderr, "too few hitpoints");

    vector<thread> ths;
    int cnt = hits.size() / UPDATE_THREADS;
    memset(cntHitUpdated, 0, sizeof(cntHitUpdated));
    for (int i = 0; i < UPDATE_THREADS; ++i) {
        auto ed = i == UPDATE_THREADS - 1 ? hits.end() : hits.begin() + (i + 1) * cnt;
        ths.push_back(thread(__updateHit, hits.begin() + i * cnt, ed, i));
    }
    int sum;
    int last = -1;
    do {
        sum = 0;
        for (int i = 0; i < UPDATE_THREADS; ++i)
            sum += cntHitUpdated[i];
        double rate = sum / (double)hits.size();
        int now = rate * 100;
        if (now != last) {
            last = now;
            printf("updating hitpoints... %d\n", now);
        }
    }
    while (sum != hits.size());
    for (auto &th : ths)
        th.join();

    delete ptntree;
    ptntree = NULL;
}

void putImage()
{
    static int id = 0;
    char str[100];
    sprintf(str, "ppm%d.png", id++);
    unsigned char *buf = (unsigned char *)malloc(WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    double *dbuf = (double *)malloc(8 * WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    memset(dbuf, 0, 8 * WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    /*static double Omega_Pixel = 4 * integral(0, camera.vx.len() / 2, [](double x) {
        static double h = (camera.o - camera.focus).len2();
        double t = 1 / sqrt(x * x + h);
        return t * atan(camera.vy.len() / 2 * t);
    });
    printf("%g\n", Omega_Pixel);*/
    double C = 1 / M_PI / Nemit / camera.vx.len() / camera.vy.len();
    for (auto &ht : hits) {
        int pos = ht.x * WINDOW_WIDTH + ht.y;
        pos *= 3;
        for (int i = 0; i < 3; ++i) {
            double tmp = ht.direct.value[i];
            if (Nemit != 0)
                tmp = (tmp + ht.tau.value[i] / sqr(ht.radius) * C);
            //if (tmp > EPS && i == 0)
            //    tmp = 100;
            //tmp *= Omega_Pixel;
            if (tmp > 255)
                tmp = 255;
            tmp *= ht.wgt.value[i];
            dbuf[pos + i] += tmp;
        }
    }
    for (int i = 0; i < WINDOW_HEIGHT * WINDOW_WIDTH * 3; ++i)
        buf[i] = dbuf[i];
    lodepng_encode24_file(str, buf, WINDOW_WIDTH, WINDOW_HEIGHT);
    free(buf);
    free(dbuf);
    printf("generated image %s\n", str);
}

void readInput()
{
    readModel("list.txt");
    ifstream fin("camera.txt");
    if (!fin.good())
        cerr << "camera.txt reading error" << endl;
    fin >> camera.focus >> camera.o >> camera.vx >> camera.vy;
    fin.close();
    camera.normal = camera.o - camera.focus;
    camera.normal /= camera.normal.len();
}

void update()
{
    static int cnt = 0;
    readfile();
    updateHitpoints();
    if (++cnt % ITERATION_PER_OUTPUT == 0)
        putImage();
}

int main(int argc, char *argv[])
{
    if (argc != 2) {
        printf("Usage: %s <directory>\n", argv[0]);
        return 1;
    }

    pid_t pid;
    if (pipe(fd2) < 0)
        cerr << "pipe error" << endl;
GO_FORK:
    if (totforks < MAXFORKS) {
        if ((pid = fork()) < 0) {
            cerr << "fork error" << endl;
            exit(1);
        }
        ++totforks;
    }
    if (pid > 0) { // parent process
        if (totforks < MAXFORKS)
            goto GO_FORK;
        close(fd2[1]);

        //puts("runnig raytracing...");
        //system((string("./raytracing ") + argv[1]).c_str());
        chdir(argv[1]);
        readInput();

        filelist.clear();
        thread th(listenTerminal);
        th.detach();

        // need to run RayTracing first
        loadHitpoints();
        putImage();

        while (!endflag)
            if (filelist.empty())
                this_thread::sleep_for(chrono::seconds(1));
            else
                update();
        while (!filelist.empty())
            update();
    }
    else {   // child process
        close(fd2[0]);
        if (fd2[1] != STDOUT_FILENO) {
            if (dup2(fd2[1], STDOUT_FILENO) != STDOUT_FILENO)
                cerr << "dup2 error to stdout" << endl;
            close(fd2[1]);
        }
        char str[10];
        sprintf(str, "%d", totforks);
        if (execl("./photontracing", "photontracing", argv[1], str, (char *)0) < 0)
            cerr << "execl error" << endl;
    }
    return 0;
}
