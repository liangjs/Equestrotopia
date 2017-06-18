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
char cwd[1000];
void make_input(int id)
{
    chdir(path.c_str());
    FILE *file = fopen("list.txt", "w");
    fprintf(file, "2\nwall.obj\nmeshs.%04d.obj\n0\n", id);
    fclose(file);
    chdir(cwd);
}
int main(int argc, char *argv[])
{
    if (argc != 2) {
        puts("Usage: dist <directory>");
        return 1;
    }
    getcwd(cwd, sizeof(cwd));
    path = argv[1];
    int L, R;
    scanf("%d%d", &L, &R);
    for (int i = L; i <= R; ++i) {
        printf("running meshs %d\n", i);

        make_input(i);

        system(("passes " + path).c_str());

        //system(("photontracing " + path).c_str());

        system(("updation " + path + " meshs" + to_string(i)).c_str());
    }
    return 0;
}
