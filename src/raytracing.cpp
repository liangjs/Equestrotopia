#include <unistd.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

#include "raytracing.h"
#include "input.h"
#include "kdtree.h"
#include "light.h"

using namespace Equestria;
using namespace std;

const double antiAliasMatrix[SAMPLE_RATE][SAMPLE_RATE] = {{0.0625, 0.125, 0.0625}, {0.125, 0.25, 0.125}, {0.0625, 0.125, 0.0625}};

polyKDTree *polytree;
Camera camera;
vector<Light> lights;

void readInput(const string &path)
{
    //cout << "reading input ..." << endl;
    chdir(path.c_str());
    readModel("list.txt");
    ifstream fin("camera.txt");
    fin >> camera.focus >> camera.o >> camera.vx >> camera.vy;
    fin.close();
    fin.open("light.txt");
    int nlight;
    fin >> nlight;
    lights.resize(nlight);
    double powersum = 0;
    for (int i = 0; i < nlight; ++i) {
        fin >> lights[i].pos >> lights[i].color >> lights[i].power;
        powersum += lights[i].power;
    }
    for (int i = 0; i < nlight; ++i)
        lights[i].power /= powersum;
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

void RayTracing(const Ray &ray, double n1, int pixel_x, int pixel_y, const Point &wgt, int deep = 0)
{
    static const int maxdeep = 3;
    if (wgt.len() < 1e-4)
        return;
    Polygon *p;
    double t = intersect(ray, polytree, p);
    if (t != INF) {
        Material::MTL &mtl = material[p->label].mtl;
        Point pos = ray.bgn + t * ray.vec;
        double R0 = sqr((n1 - mtl.Ni) / (n1 + mtl.Ni));
        double theta = dotsProduct(p->normvf, -ray.vec);
        double R = R0 + (1 - R0) * pow(1 - cos(theta), 5);
        if (deep == maxdeep || dcmp(mtl.Ns, 100) < 0) { /* not specular */
            Hitpoint hit;
            hit.position = pos;
            hit.normv = p->normvf;
            hit.raydir = -ray.vec;
            hit.material = p->label;
            hit.x = pixel_x;
            hit.y = pixel_y;
            hit.wgt = wgt * R;
            hit.direct = Point(0, 0, 0);
            for (auto &light : lights) {
                Polygon *tmp;
                double tt = intersect(Ray(light.pos, pos - light.pos), polytree, tmp);
                double vlen = (pos - light.pos).len();
                if (dcmp(tt, vlen) >= 0)
                    hit.direct += 0;
            }
            hits.push_back(hit);
        }
        else { /* specular */
            Point addition = p->normvf * EPS;
            Point dir2 = ray.vec - 2 * dotsProduct(ray.vec, p->normvf) * p->normvf;
            RayTracing(Ray(pos + addition, dir2), n1, pixel_x, pixel_y, elemMult(wgt * R, mtl.Ks), deep + 1);
            if (fabs(mtl.Tr) > EPS) { /* transparent */
                double tmp = 1 - sqr(n1) / sqr(mtl.Ni) * (1 - sqr(dotsProduct(ray.vec, p->normvf)));
                if (tmp > 0) {
                    tmp += n1 / mtl.Ni * dotsProduct(ray.vec, p->normvf);
                    dir2 = ray.vec * n1 / mtl.Ni - p->normvf * tmp;
                    RayTracing(Ray(pos + addition, dir2), mtl.Ni, pixel_x, pixel_y, elemMult(wgt * (1 - R) * mtl.Tr, mtl.Tf), deep + 1);
                }
            }
        }
    }
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
                    RayTracing(Ray(center2, center2 - camera.focus), 1, x, y, identity * antiAliasMatrix[i][j]);
                }
        }

    cout << hits.size() << endl;
    for (auto &i : hits) {
        cout << i.position << endl;
        cout << i.normv << endl;
        cout << i.material << endl;
        cout << i.x << ' ' << i.y << endl;
        cout << i.wgt << endl;
        cout << i.direct << endl;
    }
    return 0;
}
