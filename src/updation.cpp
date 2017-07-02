#include <iostream>
#include <algorithm>
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
int Nemit = 0;
char cwd[1000];
Camera camera;

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
    /*void *buf = malloc(148 * n), *ptr = buf;
    if (fread(buf, 148, n, file) != n) {
        fprintf(stderr, "hitpoints error!!!\n");
        exit(1);
    }
    fclose(file);*/
    double R0 = INITRADIUS_RATIO * sqrt(camera.vx.len() * camera.vy.len());
    for (int i = 0; i < n; ++i) {
        Hitpoint cur;
        cur.position.Read(file);
        cur.normv.Read(file);
        cur.raydir.Read(file);
        fread(&cur.mtl_label, 4, 1, file);
        fread(&cur.u, 8, 1, file);
        fread(&cur.v, 8, 1, file);
        fread(&cur.x, 4, 1, file);
        fread(&cur.y, 4, 1, file);
        /*memcpy(&cur.mtl_label, ptr, 4), ptr = (char *)ptr + 4;
        memcpy(&cur.u, ptr, 8), ptr = (char *)ptr + 8;
        memcpy(&cur.v, ptr, 8), ptr = (char *)ptr + 8;
        memcpy(&cur.x, ptr, 4), ptr = (char *)ptr + 4;
        memcpy(&cur.y, ptr, 4), ptr = (char *)ptr + 4;*/
        cur.wgt.Read(file);
        cur.direct.Read(file);
        cur.radius = R0;
        hits.push_back(cur);
    }
    //free(buf);
    fclose(file);
    random_shuffle(hits.begin(), hits.end());
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
    clearptnlist();
    const char *curfile = filelist.front().c_str();
    FILE *file = fopen(curfile, "rb");
    if (file == NULL)
        exit(1);
    Nemit += PHOTONSPER;
    int N;
    fread(&N, 4, 1, file);
    for (int i = 0; i < N; ++i) {
        Photon *curptn = new Photon;
        curptn->light.bgn.Read(file);
        curptn->light.vec.Read(file);
        curptn->rgb.Read(file);
        fread(&curptn->mtl, 4, 1, file);
        ptns.push_back(curptn);
    }
    fclose(file);
    //remove(curfile);
    filelist.pop_front();
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
            for (auto &ptn : nearptns)
                if (dotsProduct(ptn->light.vec, ht.normv) > EPS) {
                    Point brdf_cos = mt.BRDF_cos(ptn->light.vec, ht.raydir, ht.normv, ht.u, ht.v);
                    ht.tau += elemMult(ptn->rgb, brdf_cos);
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
            this_thread::sleep_for(chrono::seconds(1));
        }
    }
    while (sum != hits.size());
    for (auto &th : ths)
        th.join();

    delete ptntree;
    ptntree = NULL;
}

void putImage(const char *fname = NULL)
{
    static int id = 0;
    char str[100];
    if (fname == NULL)
        sprintf(str, "ppm%d.png", id++);
    else
        sprintf(str, "%s.png", fname);
    unsigned char *buf = (unsigned char *)malloc(WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    double *dbuf = (double *)malloc(8 * WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    memset(dbuf, 0, 8 * WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    /*static double Omega_Pixel = 4 * integral(0, camera.vx.len() / 2, [](double x) {
        static double h = (camera.o - camera.focus).len2();
        double t = 1 / sqrt(x * x + h);
        return t * atan(camera.vy.len() / 2 * t);
    });
    printf("%g\n", Omega_Pixel);*/
    double C = PHOTON_CONTRIB / M_PI / Nemit;
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
            tmp *= ht.wgt.value[i];
            dbuf[pos + i] += tmp;
        }
    }
    double mxV = 0;
    for (int i = 0; i < WINDOW_HEIGHT * WINDOW_WIDTH * 3; ++i) {
        buf[i] = min(255.0, dbuf[i]);
        mxV = max(mxV, dbuf[i]);
    }
    //printf("max value = %g\n", mxV);
    lodepng_encode24_file(str, buf, WINDOW_WIDTH, WINDOW_HEIGHT);
    free(buf);
    free(dbuf);
    printf("generated image %s\n", str);
}

void readInput()
{
    cout << "reading model..." << endl;
    readModel("list.txt");
    polygon.clear();
    ifstream fin("camera.txt");
    if (!fin.good())
        cerr << "camera.txt reading error" << endl;
    fin >> camera.focus >> camera.o >> camera.vx >> camera.vy;
    fin.close();
    camera.normal = camera.o - camera.focus;
    camera.normal /= camera.normal.len();
    fin.open("photon_contrib.txt");
    fin >> PHOTON_CONTRIB;
    fin.close();
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
    if (argc != 3) {
        printf("Usage: %s <directory> <image>\n", argv[0]);
        return 1;
    }

    //puts("runnig raytracing...");
    //system((string("./raytracing ") + argv[1]).c_str());
    chdir(argv[1]);
    readInput();

    // need to run RayTracing first
    loadHitpoints();
    putImage();

    for (int i = 0; i < MAXITERATION; ++i)
        filelist.push_back("PhotonMap" + to_string(i) + ".map");

    while (!filelist.empty())
        update();

    putImage(argv[2]);
    return 0;
}
