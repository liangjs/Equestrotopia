#include <unistd.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <fstream>

#include <GL/glut.h>

#include "photontracing.h"
#include "input.h"
#include "kdtree.h"
#include "light.h"

using namespace Equestria;
using namespace std;

polyKDTree* polytree;
Camera camera;
vector<Light> lights;

void readInput(const string& path) {
    chdir(path.c_str());
    readModel("list.txt");
    chdir("..");
}

void placeCameraAndLight() {
    double bdmin[3], bdmax[3];
    for (int i = 0; i < 3; ++i)
        bdmin[i] = INF, bdmax[i] = -INF;
    for (auto& i : polygon)
        for (int j = 0; j < 3; ++j) {
            bdmin[j] = min(bdmin[j], i->bdmin[j]);
            bdmax[j] = max(bdmax[j], i->bdmax[j]);
        }

    camera.o.value[0] = (bdmin[0] + bdmax[0]) / 2;
    camera.o.value[1] = (bdmin[1] + bdmax[1]) / 2;
    camera.o.value[2] = 2 * bdmax[2] - bdmin[2];
    camera.vx = Point(0, (bdmax[1] - bdmin[1]) / 10, 0);
    camera.vy = Point((bdmax[0] - bdmin[0]) / 10, 0, 0);
    double lvx = camera.vx.len(), lvy = camera.vy.len();
    if (lvx * WINDOW_WIDTH < lvy * WINDOW_HEIGHT / WINDOW_WIDTH)
        camera.vx *= lvy / lvx * WINDOW_HEIGHT / WINDOW_WIDTH;
    else
        camera.vy *= lvx / lvy * WINDOW_WIDTH / WINDOW_HEIGHT;
    camera.focus = camera.o;
    camera.focus.z += .11 * (bdmax[2] - bdmin[2]);
    camera.normal = camera.o - camera.focus;
    camera.normal /= camera.normal.len();

    Light l0;
    l0.pos = camera.focus;
    l0.pos.y -= (bdmax[1] - bdmin[1]) * 0.1;
    l0.color = Point(1, 1, 1);
    l0.power = 1;
    lights.push_back(l0);
}

void build_polyKDTree() {
    polytree = new polyKDTree(polygon.begin(), polygon.end());
}

Point buf[WINDOW_HEIGHT][WINDOW_WIDTH];

void ejectphoton(const Light& light, double x, double y, double z) {
    double n1 = 1;
    Photon ptn(Ray(light.pos, Point(x, y, z)), light.color);
    for (int bounce = 0; bounce < MAXBOUNCES; ++bounce) {
        Polygon* p;
        double t = intersect(ptn.light, polytree, p);
        if (t == INF) return;
        // need to consider the color change
        cout << ptn << endl;
        Material::MTL& mtl = material[p->label].mtl;
        Point pos = ptn.light.bgn + t * ptn.light.vec;
        Point N = p->getNormal(pos);
        Point addition = N * EPS;
        Point dir2 = ptn.light.vec - 2 * dotsProduct(ptn.light.vec, N) * N;
        double R0 = sqr((n1 - mtl.Ni) / (n1 + mtl.Ni));
        double R = R0 + (1 - R0) * pow(1 - dotsProduct(N, -ptn.light.vec), 5);
        if ((rand() % 10001) / 10000.0 < R) { // reflected
            ptn.light.bgn = pos + addition;
            ptn.light.vec = dir2;
        } else { // refracted or absorbed
            if (mtl.Tr > EPS) {
                double Prefraction = mtl.Tr;
                double Pabsorption = 1 - mtl.Tr;
                double tmp = 1 - sqr(n1) / sqr(mtl.Ni) * (1 - sqr(dotsProduct(ptn.light.vec, N)));
                if (tmp <= 0 || (rand() % 10001) / 10000.0 < Pabsorption)
                    return; // total reflection or absorbed
                tmp = sqrt(tmp) + n1 / mtl.Ni * dotsProduct(ptn.light.vec, N);
                dir2 = ptn.light.vec * n1 / mtl.Ni - N * tmp;
                ptn.light.bgn = pos + addition;
                ptn.light.vec = dir2;
            } else return; // absorbed
        }
    }
}

unsigned char pixels[WINDOW_HEIGHT * WINDOW_WIDTH * 3];

void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(WINDOW_WIDTH, WINDOW_HEIGHT, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    glFlush();
    glutSwapBuffers();
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "usage: photontracing model_directory" << endl;
        return -1;
    }

    srand(time(0));
    readInput(argv[1]);
    placeCameraAndLight();
    build_polyKDTree();

    for (int iteration = 0; iteration < MAXITERATION; ++iteration) {
        for (auto& curlight : lights) {
            for (int i = 0; i < PHOTONSPER * curlight.power; ++i) {
                double theta = (rand() % 10000) / 5000.0 * M_PI;
                double z = (rand() % 10001) / 5000.0 - 1;
                double x = sqrt(1 - sqr(z)) * cos(theta),
                       y = sqrt(1 - sqr(z)) * sin(theta);
                ejectphoton(curlight, x, y, z);
            }
        }
        printf("Iteration %d Done\n", iteration);
    }

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    glutCreateWindow("PhotonTracing");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMainLoop();
    return 0;
}
