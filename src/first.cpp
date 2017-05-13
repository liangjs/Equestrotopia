#include <unistd.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "first.h"
#include "input.h"
#include "kdtree.h"
#include "light.h"

using namespace Equestria;
using namespace std;

const double antiAliasMatrix[SAMPLE_RATE][SAMPLE_RATE] = {{0.0625, 0.125, 0.0625}, {0.125, 0.25, 0.125}, {0.0625, 0.125, 0.0625}};

polyKDTree *polytree;
Camera camera;

void readInput(const string &path)
{
    //cout << "reading input ..." << endl;
    chdir(path.c_str());
    readModel("list.txt");
    ifstream fin("camera.txt");
    fin >> camera.focus >> camera.o >> camera.vx >> camera.vy;
    fin.close();
    camera.normal = camera.o - camera.focus;
    camera.normal /= camera.normal.len();
    chdir("..");
}

void build_polyKDTree()
{
    //cout << "building polygon kd-tree ..." << endl;
    polytree = new polyKDTree(polygon.begin(), polygon.end());
}

std::vector<Hitpoint> hits;

void RayTracing(const Ray &ray, int pixel_x, int pixel_y, const Point &wgt, int deep = 0)
{
    static const int maxdeep = 3;
}

int main(int argc, char *argv[])
{
    if (argc == 1 || argc > 2) {
        cout << "usage: Equestrotopia model_directory" << endl;
        return -1;
    }

    readInput(argv[1]);
    build_polyKDTree();

    Point dx = camera.vx / WINDOW_HEIGHT;
    Point dy = camera.vy / WINDOW_WIDTH;
    Point _dx = dx / SAMPLE_RATE, _dy = dy / SAMPLE_RATE;
    Point identity(1, 1, 1);
    for (int x = 0; x < WINDOW_HEIGHT; ++x)
        for (int y = 0; y < WINDOW_WIDTH; ++y) {
            Point center = camera.o + (x - WINDOW_HEIGHT / 2.0 + 0.5) * dx + (y - WINDOW_WIDTH / 2.0 + 0.5) * dy;
            for (int i = 0; i < SAMPLE_RATE; ++i)
                for (int j = 0; j < SAMPLE_RATE; ++j) {
                    Point center2 = center + (i - SAMPLE_RATE / 2.0 + 0.5) * _dx + (j - SAMPLE_RATE / 2.0 + 0.5) * _dy;
                    RayTracing(Ray(center2, camera.normal), x, y, identity * antiAliasMatrix[i][j]);
                }
        }

    cout << hits.size() << endl;
    for (auto &i : hits) {
        cout << i.position << endl;
        cout << i.normv << endl;
        cout << i.material << endl;
        cout << i.x << ' ' << i.y << endl;
        cout << i.wgt << endl;
    }
    return 0;
}
