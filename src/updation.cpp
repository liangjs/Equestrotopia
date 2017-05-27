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
#include "updation.h"

using namespace std;
using namespace Equestria;

list<FILE*> filelist;
vector<Hitpoint> hits;
int totforks = 0;
int fd1[2], fd2[2];
mutex filelistlock;
bool endflag = false;
streambuf* coutBuf = cout.rdbuf();
streambuf* cinBuf = cin.rdbuf();

void sig_pipe(int signo) {
    printf("SIGPIPE caught\n");
    exit(1);
}

void listenTerminal() {
    char line[MAXNLINE];
    int n;
    while ((n = read(fd2[0], line, MAXNLINE)) >= 0) {
        if (n == 0) {
            printf("child closed pipe\n");
            break;
        }
        char str[100];
        sscanf(line, "%s", str);
        filelistlock.lock();
        filelist.push_back(fopen(str, "r"));
        filelistlock.unlock();
    }
    if (ferror(stdin))
        cerr << "fgets error on stdin" << endl;
    endflag = true;
}

void loadHitpoints() {
    FILE* input = fopen("hitpoints.txt", "r");
    ifstream iff("hitpoints.txt");
    cin.rdbuf(iff.rdbuf());
    int n;
    fscanf(input, "%d", &n);
    for (int i = 0; i < n; ++i) {
        Hitpoint cur;
        cin >> cur.position >> cur.normv
            >> cur.material
            >> cur.x >> cur.y
            >> cur.wgt
            >> cur.direct;
        hits.push_back(cur);
    }
    iff.close();
    cin.rdbuf(cinBuf);
    fclose(input);
}

void updateHitpoints() {

}

int main(int argc, char* argv[]) {
    pid_t pid;

    if (signal(SIGPIPE, sig_pipe) == SIG_ERR)
        cerr << "signal error" << endl;
    if (pipe(fd1) < 0 || pipe(fd2) < 0)
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
        close(fd1[0]);
        close(fd1[1]);
        close(fd2[1]);

        filelist.clear();
        thread th(listenTerminal);
        th.detach();

        // need to run RayTracing first
        loadHitpoints();

        while (!endflag)
            if (!filelist.empty())
                updateHitpoints();

        exit(0);
    } else { // child process
        close(fd1[0]);
        close(fd1[1]);
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
