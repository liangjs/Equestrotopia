#include "dist.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>
#include <unistd.h>
#include <thread>
#include <mutex>
#include <condition_variable>
using namespace std;
string path;
//char cwd[1000];
void emit_photons()
{
    mutex runMutex;
    condition_variable runCV;
    int resource = RUNTHREADS;
    auto run = [&](int id) {
        {
            unique_lock<mutex> runLock(runMutex);
            if (resource > 0)
                --resource;
            else {
                runCV.wait(runLock, [&] {return resource > 0;});
                --resource;
            }
        }
        system(("./photontracing " + path + " " + to_string(id)).c_str());
        {
            unique_lock<mutex> runLock(runMutex);
            ++resource;
            runCV.notify_one();
        }
    };
    vector<thread> ths;
    for (int i = 0; i < PHOTON_FORKS; ++i)
        ths.push_back(thread(run, i));
    for (auto &th : ths)
        th.join();
}
void make_input(int id)
{
    FILE *file = fopen("list.txt", "w");
    fprintf(file, "2\nwall.obj\nmeshs.%04d.obj\n0\n", id);
    fclose(file);
}
int main(int argc, char *argv[])
{
    if (argc != 2) {
        puts("Usage: dist <directory>");
        return 1;
    }
    //getcwd(cwd, sizeof(cwd));
    path = argv[1];
    for (int i = 0; i < NMESHS; ++i) {
        printf("running (%d)\n", i);

        make_input(i);

        system(("./raytracing " + path).c_str());

        emit_photons();

        system(("./updation " + path).c_str());
    }
    return 0;
}
