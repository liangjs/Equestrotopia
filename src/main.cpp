#include <unistd.h>
#include <string>
#include <iostream>
#include <algorithm>

#include "main.h"
#include "input.h"
#include "kdtree.h"

using namespace Equestria;
using namespace std;
polyKDTree *polytree;

void readInput(const string &path)
{
    cout << "reading input ..." << endl;

    string spath = path;
    spath.erase(spath.rfind('/'));
    chdir(spath.c_str());

    chdir("data");
    readModel("list.txt");
    chdir("..");
}

void build_polyKDTree()
{
    cout << "building polygon kd-tree ..." << endl;

    polytree = new polyKDTree(polygon.begin(), polygon.end());
}

int main(int argc, char *argv[])
{
    readInput(argv[0]);
    build_polyKDTree();
    return 0;
}
