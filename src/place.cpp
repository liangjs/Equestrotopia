#include <unistd.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "raytracing.h"
#include "input.h"
#include "light.h"

using namespace Equestria;
using namespace std;

polyKDTree *polytree;
Camera camera;
vector<Light> lights;

void readInput()
{
    readModel("list.txt");
    rotateModel("prerotate.txt");
}

void placeCameraLight()
{
    double bdmin[3], bdmax[3];
    for (int i = 0; i < 3; ++i)
        bdmin[i] = INF, bdmax[i] = -INF;
    auto wallit = mtlIndex.find("wall");
    int mtlWall;
    if (wallit == mtlIndex.end())
        mtlWall = -1;
    else
        mtlWall = wallit->second;
    for (auto &i : polygon)
        if (i->label != mtlWall)
            for (int j = 0; j < 3; ++j) {
                bdmin[j] = min(bdmin[j], i->bdmin[j]);
                bdmax[j] = max(bdmax[j], i->bdmax[j]);
            }

    Point center((bdmin[0] + bdmax[0]) / 2, (bdmin[1] + bdmax[1]) / 2, (bdmin[2] + bdmax[2]) / 2);

    camera.o = center;
    camera.o.y = bdmin[1] * 0.9 + bdmin[1] * 0.1;
    camera.vx = Point(0, 0, -(bdmax[2] - bdmin[2]) / 2);
    camera.vy = Point((bdmax[0] - bdmin[0]) / 2, 0, 0);
    double lvx = camera.vx.len(), lvy = camera.vy.len();
    if (lvx * WINDOW_WIDTH < lvy * WINDOW_HEIGHT)
        camera.vx *= lvy / lvx * WINDOW_HEIGHT / WINDOW_WIDTH;
    else
        camera.vy *= lvx / lvy * WINDOW_WIDTH / WINDOW_HEIGHT;
    camera.focus = center;
    camera.focus.y -= bdmax[1] - bdmin[1];
    camera.normal = camera.o - camera.focus;
    camera.normal /= camera.normal.len();

    Light l0, l1;
    l0.pos = camera.o;
    l0.pos.z += (bdmax[2] - bdmin[2]) / 2;
    l0.color = Point(1, 1, 1);
    l0.power = 0.5;
    lights.push_back(l0);

    l1.pos = center;
    l1.pos.z += (bdmax[2] - bdmin[2]) / 2;
    l1.color = Point(1, 1, 1);
    l1.power = 0.5;
    lights.push_back(l1);
}

int main(int argc, char *argv[])
{
    if (argc == 1 || argc > 2) {
        cout << "usage: place model_directory" << endl;
        return -1;
    }
    chdir(argv[1]);

    readInput();
    placeCameraLight();

    ofstream fout("camera.txt");
    fout << camera.focus << endl << camera.o << endl << camera.vx << ' ' << camera.vy << endl;
    fout.close();

    fout.open("light.txt");
    fout << lights.size() << endl;
    for (auto &i : lights)
        fout << i.pos << ' ' << i.color << ' ' << i.power << endl;
    fout.close();

    return 0;
}
