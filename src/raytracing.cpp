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

Camera camera;
vector<Light> lights;
std::vector<Hitpoint> hits;

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

void RayTracing(const Ray &ray, double n1, int pixel_x, int pixel_y, const Point &wgt, int deep = 0)
{
    static const int maxdeep = 3;
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
            RayTracing(r2, n1, pixel_x, pixel_y, elemMult(wgt, brdf_cos), deep + 1);
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
        hits.push_back(hit);
    }
}

void Run()
{
    Point dx = camera.vx / WINDOW_HEIGHT;
    Point dy = camera.vy / WINDOW_WIDTH;
    Point _dx = dx / SAMPLE_RATE, _dy = dy / SAMPLE_RATE;
    Point identity(1, 1, 1);
    for (int x = 0; x < WINDOW_HEIGHT; ++x) {
        //printf("progress: %g\n", (float)x / WINDOW_WIDTH);
        for (int y = 0; y < WINDOW_WIDTH; ++y) {
            Point center = camera.o + (x - WINDOW_HEIGHT / 2.0 + 0.5) * dx + (y - WINDOW_WIDTH / 2.0 + 0.5) * dy;
            for (int i = 0; i < SAMPLE_RATE; ++i)
                for (int j = 0; j < SAMPLE_RATE; ++j) {
                    Point center2 = center + (i - SAMPLE_RATE / 2.0 + 0.5) * _dx + (j - SAMPLE_RATE / 2.0 + 0.5) * _dy;
                    RayTracing(Ray(center2, center2 - camera.focus), 1, x, y, identity * antiAliasMatrix[i][j]);
                }
        }
    }
}

#include "lodepng.h"

void output()
{
    /*FILE *fppm = fopen("test.ppm", "wb");
    fprintf(fppm, "P6\n%d %d\n255\n", WINDOW_WIDTH, WINDOW_HEIGHT);
    unsigned char *buf = (unsigned char *)malloc(WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    double *dbuf = (double *)malloc(8*WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    memset(dbuf, 0, 8 * WINDOW_WIDTH * WINDOW_HEIGHT * 3);
    for (auto &ht : hits) {
        int pos = ht.x * WINDOW_WIDTH + ht.y;
        pos *= 3;
        for (int i = 0; i < 3; ++i) {
            double tmp = ht.direct.value[i];
            tmp *= 255;
            if (tmp > 255)
                tmp = 255;
            tmp *= ht.wgt.value[i];
            dbuf[pos + i] += tmp;
        }
    }
    for (int i = 0; i < WINDOW_WIDTH * WINDOW_HEIGHT * 3; ++i)
        buf[i] = dbuf[i];
    lodepng_encode24_file("test.png", buf, WINDOW_WIDTH, WINDOW_HEIGHT);
    fwrite(buf, 1, WINDOW_HEIGHT * WINDOW_WIDTH * 3, fppm);
    free(buf);
    free(dbuf);
    fclose(fppm);*/

    FILE *file = fopen("hitpoints.txt", "wb");
    int sz = hits.size();
    fwrite(&sz, 4, 1, file);
    for (auto &i : hits) {
        i.position.Print(file);
        i.normv.Print(file);
        i.raydir.Print(file);
        fwrite(&i.mtl_label, 4, 1, file);
        fwrite(&i.u, 8, 1, file);
        fwrite(&i.v, 8, 1, file);
        fwrite(&i.x, 4, 1, file);
        fwrite(&i.y, 4, 1, file);
        i.wgt.Print(file);
        i.direct.Print(file);
    }
    fclose(file);
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

    Run();

    output();
    return 0;
}
