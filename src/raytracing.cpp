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
std::vector<Hitpoint> hits;

void readInput(const string &path)
{
    //cout << "reading input ..." << endl;
    chdir(path.c_str());
    readModel("list.txt");
    rotateModel("prerotate.txt");
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

void RayTracing(const Ray &ray, double n1, int pixel_x, int pixel_y, const Point &wgt, int deep = 0)
{
    static const int maxdeep = 3;
    if (wgt.len() < 1e-4)
        return;
    Polygon *p;
    double t = intersect(ray, polytree, p);
    if (t == INF)
        return;
    Material::MTL &mtl = material[p->label].mtl;
    Point pos = ray.bgn + t * ray.vec;
    Point N = p->getNormal(pos);
    double Ni = mtl.Ni;
    if (dotsProduct(N, -ray.vec) < 0) { /* go out of the material */
        Ni = 1;
        N = -N;
    }
    bool end = true;
    if (deep < maxdeep) {
        if (dcmp(mtl.Ns, 100) > 0) {/* specular */
            RayTracing(reflect(pos, N, ray.vec), n1, pixel_x, pixel_y, elemMult(wgt, mtl.Ks), deep + 1);
            end = false;
        }
        if (mtl.Tr > EPS) /* transparent */
            try {
                RayTracing(refract(pos, n1, Ni, N, ray.vec), Ni, pixel_x, pixel_y, elemMult(mtl.Tf, wgt * mtl.Tr), deep + 1);
                end = false;
            }
            catch (std::logic_error) {}
    }
    if (end) {
        Hitpoint hit;
        hit.position = pos;
        hit.normv = N;
        hit.raydir = -ray.vec;
        hit.material = p->label;
        hit.x = pixel_x;
        hit.y = pixel_y;
        hit.wgt = wgt;
        hit.direct = Point(0, 0, 0);
        for (auto &light : lights) {
            Polygon *tmp;
            Point L = light.pos - pos;
            if (dotsProduct(L, N) <= EPS) // inside the object, no direct light
                continue;
            double tt = intersect(Ray(light.pos, -L), polytree, tmp);
            double vlen = L.len();
            if (dcmp(tt, vlen) >= 0) {
                L /= vlen;
                Point c = elemMult(light.color, material[p->label].BRDF(L, -ray.vec, N));
                hit.direct += c * light.power;
            }
        }
        hits.push_back(hit);
    }
}

int main(int argc, char *argv[])
{
    if (argc == 1 || argc > 2) {
        cout << "usage: raytracing model_directory" << endl;
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

    ofstream of("hitpoints.txt");
    of << hits.size() << endl;
    for (auto &i : hits) {
        of << i.position << endl;
        of << i.normv << endl;
        of << i.raydir << endl;
        of << i.material << endl;
        of << i.x << ' ' << i.y << endl;
        of << i.wgt << endl;
        of << i.direct << endl;
    }
    of.close();
    return 0;
}
