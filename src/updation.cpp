#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <cstring>
#include <list>
#include <vector>
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <thread>
#include <mutex>

#include "light.h"
#include "input.h"
#include "kdtree.h"
#include "updation.h"

#include "lodepng.h"

using namespace std;
using namespace Equestria;

ptnKDTree* ptntree = NULL;
list<string> filelist;
vector<Hitpoint> hits;
vector<Photon*> ptns;
vector<Photon*> nearptns;
int totforks = 0;
int fd1[2], fd2[2];
mutex filelistlock;
bool endflag = false;
streambuf* coutBuf = cout.rdbuf();
streambuf* cinBuf = cin.rdbuf();
int Nemit = 0;

void pushFile(const char* str) {
    filelistlock.lock();
    filelist.push_back(str);
    filelistlock.unlock();
}

void listenTerminal() {
    char line[MAXNLINE];
    int n;
    FILE *file = fdopen(fd2[0], "r");
    while (fscanf(file, "%s", line) != EOF)
        pushFile(line);
    endflag = true;
}

void loadHitpoints() {
    ifstream iff("hitpoints.txt");
    cin.rdbuf(iff.rdbuf());
    int n;
    cin >> n;
    for (int i = 0; i < n; ++i) {
        Hitpoint cur;
        cin >> cur.position >> cur.normv
            >> cur.material
            >> cur.raydir
            >> cur.x >> cur.y
            >> cur.wgt
            >> cur.direct;
        cur.radius = INITRADIUS;
        hits.push_back(cur);
    }
    iff.close();
    cin.rdbuf(cinBuf);
}

void clearptnlist() {
    while (!ptns.empty()) {
        delete ptns.back();
        ptns.pop_back();
    }
}

void readfile() {
    filelistlock.lock();
    clearptnlist();
    const char* curfile = filelist.front().c_str();
    ifstream iff(curfile);
    cin.rdbuf(iff.rdbuf());
    filelist.pop_front();
    Photon* curptn = new Photon;
    while (cin >> curptn->light.bgn
            >> curptn->light.vec
            >> curptn->rgb
            >> curptn->mtl) {
        ptns.push_back(curptn);
        curptn = new Photon;
    }
    delete curptn;
    iff.close();
    cin.rdbuf(cinBuf);
    filelistlock.unlock();
}

void updateHitpoints() {
    readfile();
    ptntree = new ptnKDTree(ptns.begin(), ptns.end());

    for (auto& ht : hits)
        if (ht.radius > EPS) {
            nearptns.clear();
            ptntree->find(ht, nearptns);
            int M = nearptns.size();
            double N = ht.ptncount;
            ht.ptncount += ALPHA * M;
            double ratio = ht.ptncount / (N + M);
            ht.radius *= sqrt(ratio);
            for (auto &ptn : nearptns)
                ht.tau += elemMult(ptn->rgb, material[ptn->mtl].BRDF(ptn->light.vec, ht.raydir, ht.normv));
            ht.tau *= ratio;
        }
    Nemit += ptns.size();

    delete ptntree;
    ptntree = NULL;
}

void putImage(int id) {
    char str[100];
    sprintf(str, "ppm%d.png", id);
    std::vector<std::uint8_t> Buffer(WINDOW_WIDTH * WINDOW_HEIGHT * 4); // RGBA
    for (auto &ht : hits) {
        int pos = ht.x * WINDOW_WIDTH + ht.y;
        for (int i = 0; i < 3; ++i) {
            double tmp = ht.tau.value[i] / Nemit;
            Buffer[pos + i] = uint8_t(255 * tmp);
        }
        Buffer[pos + 4] = 0;
    }
    std::vector<std::uint8_t> ImageBuffer;
    lodepng::encode(ImageBuffer, Buffer, WINDOW_WIDTH, WINDOW_HEIGHT);
    lodepng::save_file(ImageBuffer, str);
}

int main(int argc, char* argv[]) {
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

        filelist.clear();
        thread th(listenTerminal);
        th.detach();

        // need to run RayTracing first
        loadHitpoints();

        while (!endflag)
            if (!filelist.empty()) {
                static int cnt = 0;
                updateHitpoints();
                ++cnt;
                if (cnt % ITERATION_PER_OUTPUT)
                    putImage(cnt / ITERATION_PER_OUTPUT);
            }

        exit(0);
    } else { // child process
        close(fd2[0]);
        if (fd2[1] != STDOUT_FILENO) {
            if (dup2(fd2[1], STDOUT_FILENO) != STDOUT_FILENO)
                cerr << "dup2 error to stdout" << endl;
            close(fd2[1]);
        }
        char str[10];
        sprintf(str, "%d", totforks);
        if (execl("./photontracing", "photontracing", str, (char*)0) < 0)
            cerr << "execl error" << endl;
    }
    exit(0);
}
