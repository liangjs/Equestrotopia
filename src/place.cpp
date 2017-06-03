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

    camera.o.value[0] = (bdmin[0] + bdmax[0]) / 2;
    camera.o.value[1] = (bdmin[1] + bdmax[1]) / 2;
    camera.o.value[2] = 2 * bdmin[2] - bdmax[2];
    camera.vx = Point(0, (bdmax[1] - bdmin[1]) / 10, 0);
    camera.vy = Point((bdmax[0] - bdmin[0]) / 10, 0, 0);
    double lvx = camera.vx.len(), lvy = camera.vy.len();
    if (lvx * WINDOW_WIDTH < lvy * WINDOW_HEIGHT)
        camera.vx *= lvy / lvx * WINDOW_HEIGHT / WINDOW_WIDTH;
    else
        camera.vy *= lvx / lvy * WINDOW_WIDTH / WINDOW_HEIGHT;
    camera.focus = camera.o;
    camera.focus.z -= .11 * (bdmax[2] - bdmin[2]);

    Light l0, l1;
    /*l0.pos.value[0] = (bdmin[0] + bdmax[0]) / 2;
    l0.pos.value[1] = (bdmin[1] + bdmax[1]) / 2;
    l0.pos.value[2] = 2 * bdmax[2] - bdmin[2];
    //l0.pos.x += 0.1 * (bdmax[0] - bdmin[0]);
    //l0.pos.y += 0.1 * (bdmax[1] - bdmin[1]);
    l0.pos.z += 0.01 * (bdmax[2] - bdmin[2]);*/
    l0.pos = camera.focus;
    l0.pos.y += (bdmax[1] - bdmin[1]) * 0.25;
    l0.color = Point(1, 1, 1);
    l0.power = 0.5;
    lights.push_back(l0);

    l1.pos.x = (bdmin[0] + bdmax[0]) / 2;
    l1.pos.y = (bdmin[1] + bdmax[1]) / 2;
    l1.pos.z = (bdmax[2] + bdmin[2]) / 2;
    l1.pos.y -= bdmax[1] - bdmin[1];
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
