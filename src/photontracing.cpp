#include <unistd.h>
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <sys/types.h>
#include <unistd.h>

#include <GL/glut.h>

#include "photontracing.h"
#include "input.h"
#include "kdtree.h"
#include "light.h"

using namespace Equestria;
using namespace std;

Camera camera;
vector<Light> lights;
int curv;

void readInput(const string &path)
{
    chdir(path.c_str());
    readModel("list.txt");
    //rotateModel("prerotate.txt");
    ifstream fin("camera.txt");
    fin >> camera.focus >> camera.o >> camera.vx >> camera.vy;
    fin.close();
    camera.normal = camera.o - camera.focus;
    camera.normal /= camera.normal.len();
    fin.open("light.txt");
    int nlight;
    fin >> nlight;
    lights.resize(nlight);
    for (int i = 0; i < nlight; ++i)
        fin >> lights[i].pos >> lights[i].I >> lights[i].power;
    fin.close();
}

void build_polyKDTree()
{
    polytree = new polyKDTree(polygon.begin(), polygon.end());
}

Point getSphereRandomPoint(int flag = 0, Point *norm = NULL)
{
    // 0 : whole sphere
    // 1 : semi-sphere
    double theta = (rand() % 10000) / 5000.0 * M_PI;
    double z = rand() % 10001;
    if (flag == 0)
        (z /= 5000.0) -= 1;
    else if (flag == 1)
        z /= 10000.0;
    double x = sqrt(1 - sqr(z)) * cos(theta),
           y = sqrt(1 - sqr(z)) * sin(theta);
    if (!norm || flag == 0)
        return Point(x, y, z);
    Point tx;
    if (fabs(norm->x) > EPS)
        tx = Point(-norm->y, norm->x, 0);
    else
        tx = Point(0, -norm->z, norm->y);
    tx.normalize();
    Point ty = crossProduct(tx, *norm);
    return z * *norm + x * tx + y * ty;
}

void ejectphoton(const Light &light, const Point &dir, FILE *file)
{
    double n1 = 1;
    Point rgb(light.I);
    Ray ray(light.pos, dir);
    for (int bounce = 0; bounce < MAXBOUNCES; ++bounce) {
        Object *p;
        Point pos, N;
        double t = intersect(ray, pos, N, p);
        if (t == INF)
            return;
        if (bounce) { // exclude the first hit
            pos.Print(file);
            (-ray.vec).Print(file);
            rgb.Print(file);
            fwrite(&p->label, 4, 1, file);
        }
        Material::MTL &mtl = material[p->label].mtl;
        double Ni = mtl.Ni;
        if (dotsProduct(N, -ray.vec) < 0) { // inside object
            Ni = 1;
            N = -N;
        }
        double R0 = sqr((n1 - Ni) / (n1 + Ni));
        double R = R0 + (1 - R0) * pow(1 - dotsProduct(N, -ray.vec), 5);
        double u, v;
        p->txCoordinate(pos, u, v);
        if ((rand() % 10001) / 10000.0 < R) { // reflected
            ray.bgn = pos + N * EPS;
            Point reflectv = getSphereRandomPoint(1, &N);
            Point brdf_cos = material[p->label].BRDF_cos(-ray.vec, reflectv, N, u, v);
            rgb.multiByChannel(brdf_cos);
            ray.vec = reflectv;
        }
        else {   // refracted or absorbed
            if (mtl.Tr > EPS) {
                double Prefraction = mtl.Tr;
                double Pabsorption = 1 - mtl.Tr;
                if ((rand() % 10001) / 10000.0 < Pabsorption)
                    return; // absorbed
                try {
                    ray = refract(pos, n1, Ni, N, ray.vec);
                    rgb.multiByChannel(mtl.Tf);
                }
                catch (std::logic_error) {
                    return; // total inner reflect
                }
            }
            else
                return; // absorbed
        }
    }
}

int main(int argc, char *argv[])
{
    if (argc != 3) {
        printf("Usage: %s <directory> <id>\n", argv[0]);
        return 1;
    }
    sscanf(argv[2], "%d", &curv);

    srand(time(0)*getpid());
    readInput(argv[1]);
    build_polyKDTree();

    for (int iteration = 0; iteration < MAXITERATION; ++iteration) {
        char str[100];
        sprintf(str, "PhotonMap%d-%d.map", curv, iteration);
        FILE *file = fopen(str, "w");
        int Nemit = PHOTONSPER;
        fwrite(&Nemit, 4, 1, file);
        for (auto &curlight : lights)
            for (int i = 0; i < curlight.power * PHOTONSPER; ++i)
                ejectphoton(curlight, getSphereRandomPoint(), file);
        fclose(file);
        puts(str);
        fflush(stdout);
    }
    return 0;
}
