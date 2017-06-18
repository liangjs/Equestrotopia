#include <unistd.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <thread>
#include <mutex>

#include "raytracing.h"
#include "input.h"
#include "kdtree.h"
#include "light.h"

using namespace Equestria;
using namespace std;

Camera camera;
vector<Light> lights;
//std::vector<Hitpoint> hits;

FILE *fout;
mutex fmutex;

void readInput()
{
    //cout << "reading input ..." << endl;
    readModel("list.txt");
    //rotateModel("prerotate.txt");
    ifstream fin("camera.txt");
    if (!fin.good())
        cerr << "camera.txt reading error" << endl;
    fin >> camera.focus >> camera.o >> camera.vx >> camera.vy;
    fin.close();
    fin.open("light.txt");
    if (!fin.good())
        cerr << "light.txt reading error" << endl;
    int nlight;
    fin >> nlight;
    lights.resize(nlight);
    for (int i = 0; i < nlight; ++i)
        fin >> lights[i].pos >> lights[i].I >> lights[i].power;
    fin.close();
    camera.normal = camera.o - camera.focus;
    camera.normal /= camera.normal.len();
}

void build_polyKDTree()
{
    //cout << "building polygon kd-tree ..." << endl;
    polytree = new polyKDTree(polygon.begin(), polygon.end());
}

int nhit;

void RayTracing(const Ray &ray, double n1, int pixel_x, int pixel_y, const Point &wgt, int deep = 0)
{
    static const int maxdeep = 5;
    if (wgt.len() < 1e-4)
        return;
    Object *p;
    Point pos, N;
    double t = intersect(ray, pos, N, p);
    if (t == INF)
        return;
    int label = p->label;
    Material::MTL &mtl = material[label].mtl;
    double Ni = mtl.Ni;
    if (dotsProduct(N, -ray.vec) < 0) { /* go out of the material */
        Ni = 1;
        N = -N;
    }
    bool end = true;
    double u, v;
    p->txCoordinate(pos, u, v);
    Point Kd = mtl.getKd(u, v), Ks = mtl.getKs(u, v);
    if (deep < maxdeep) {
        if (dcmp(mtl.Ns, 100) > 0) {/* specular */
            Ray r2(reflect(pos, N, ray.vec));
            Point brdf_cos = material[p->label].BRDF_cos(r2.vec, -ray.vec, N, u, v);
            RayTracing(r2, n1, pixel_x, pixel_y, elemMult(wgt, brdf_cos) * (1 - mtl.Tr), deep + 1);
            end = false;
        }
        if (mtl.Tr > EPS) /* transparent */
            try {
                RayTracing(refract(pos, n1, Ni, N, ray.vec), Ni, pixel_x, pixel_y, elemMult(mtl.Tf, wgt * mtl.Tr), deep + 1);
                end = false;
            }
            catch (std::logic_error) {}
    }
    if (deep * 2 > maxdeep || end) {
        Hitpoint hit;
        hit.position = pos;
        hit.normv = N;
        hit.raydir = -ray.vec;
        hit.mtl_label = label;
        hit.u = u, hit.v = v;
        hit.x = pixel_x;
        hit.y = pixel_y;
        hit.wgt = wgt;
        hit.direct = Point(0, 0, 0);
        for (auto &light : lights) {
            Object *tmp;
            Point L = light.pos - pos;
            if (dotsProduct(L, N) <= EPS) // inside the object, no direct light
                continue;
            Point _pos, _N;
            double tt = intersect(Ray(light.pos, -L), _pos, _N, tmp);
            double vlen = L.len();
            if (dcmp(tt, vlen) >= 0) {
                L /= vlen;
                Point brdf_cos = material[label].BRDF_cos(L, -ray.vec, N, u, v);
                hit.direct += elemMult(light.I / sqr(vlen) * light.power, brdf_cos);
            }
        }
        fmutex.lock();
        hit.position.Print(fout);
        hit.normv.Print(fout);
        hit.raydir.Print(fout);
        fwrite(&hit.mtl_label, 4, 1, fout);
        fwrite(&hit.u, 8, 1, fout);
        fwrite(&hit.v, 8, 1, fout);
        fwrite(&hit.x, 4, 1, fout);
        fwrite(&hit.y, 4, 1, fout);
        hit.wgt.Print(fout);
        hit.direct.Print(fout);
        ++nhit;
        fmutex.unlock();
        //hits.push_back(hit);
    }
}

void Run()
{
    Point dx = camera.vx / WINDOW_HEIGHT;
    Point dy = camera.vy / WINDOW_WIDTH;
    Point _dx = dx / SAMPLE_RATE, _dy = dy / SAMPLE_RATE;
    Point identity(1, 1, 1);
    struct T {
        Ray ray;
        int x, y;
        Point wgt;
        T(const Ray &_ray, int _x, int _y, const Point &_wgt)
            : ray(_ray), x(_x), y(_y), wgt(_wgt) {}
    };
    vector<T> v[RAYTRACING_THREADS];
    for (int x = 0; x < WINDOW_HEIGHT; ++x) {
        for (int y = 0; y < WINDOW_WIDTH; ++y) {
            Point center = camera.o + (x - WINDOW_HEIGHT / 2.0 + 0.5) * dx + (y - WINDOW_WIDTH / 2.0 + 0.5) * dy;
            for (int i = 0; i < SAMPLE_RATE; ++i)
                for (int j = 0; j < SAMPLE_RATE; ++j) {
                    Point center2 = center + (i - SAMPLE_RATE / 2.0 + 0.5) * _dx + (j - SAMPLE_RATE / 2.0 + 0.5) * _dy;
                    static int cnt = 0;
                    v[cnt++].push_back(T(Ray(center2, center2 - camera.focus), x, y, identity * antiAliasMatrix[i][j]));
                    cnt %= RAYTRACING_THREADS;
                }
        }
    }

    vector<thread> ths;
    int cntRayTraced[RAYTRACING_THREADS] = {};
    auto run = [&](int id, const vector<T> &v) {
        for (auto &t : v) {
            RayTracing(t.ray, 1, t.x, t.y, t.wgt);
            ++cntRayTraced[id];
        }
    };
    for (int i = 0; i < RAYTRACING_THREADS; ++i)
        ths.push_back(thread(run, i, v[i]));
    int sum;
    int last = -1;
    const int SUM = WINDOW_HEIGHT * WINDOW_WIDTH * SAMPLE_RATE * SAMPLE_RATE;
    do {
        sum = 0;
        for (int i = 0; i < RAYTRACING_THREADS; ++i)
            sum += cntRayTraced[i];
        double rate = sum / (double)SUM;
        int now = rate * 100;
        if (now != last) {
            last = now;
            printf("raytracing... %d\n", now);
            this_thread::sleep_for(chrono::seconds(1));
        }
    }
    while (sum != SUM);
    for (auto &th : ths)
        th.join();
}

int main(int argc, char *argv[])
{
    if (argc == 1 || argc > 2) {
        cout << "usage: raytracing model_directory" << endl;
        return -1;
    }

    chdir(argv[1]);
    readInput();
    build_polyKDTree();

    fout = fopen("hitpoints.txt", "wb");
    fseek(fout, 4, SEEK_SET);

    Run();

    fseek(fout, 0, SEEK_SET);
    fwrite(&nhit, 4, 1, fout);
    fclose(fout);
    return 0;
}
