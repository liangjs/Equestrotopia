#include <GL/glut.h>
#include <unistd.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "main.h"
#include "kdtree.h"

Equestria::polyKDTree *ptree;

void display()
{
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0, 0, 1.5, 0, 0, 0, 0, 1, 0);
    {
        GLfloat light_pos[] = {0, 0, 1.1, 1};
        GLfloat light_am[] = {0.3, 0.3, 0.3, 1};
        GLfloat light_df[] = {0.9, 0.9, 0.9, 1};
        GLfloat light_sp[] = {1, 1, 1, 1};
        //glTranslatef(light_pos[0], light_pos[1], light_pos[2]);
        //glutSolidSphere(0.1, 50, 50);
        //glTranslatef(-light_pos[0], -light_pos[1], -light_pos[2]);
        glLightfv(GL_LIGHT0, GL_POSITION, light_pos);
        glLightfv(GL_LIGHT0, GL_AMBIENT, light_am);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_df);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_sp);
    }
    //glRotated(90, 0.3, 0.2, 0.7);
    static double angle = 0;
    glRotated(angle += 0.5, 0, 1, 0);
    /*{
      glBegin(GL_LINES);
      glColor3f(1, 0, 0);
      glVertex3f(0, 0, 0);
      glVertex3f(1, 0, 0);
      glColor3f(0, 1, 0);
      glVertex3f(0, 0, 0);
      glVertex3f(0, 1, 0);
      glColor3f(0, 0, 1);
      glVertex3f(0, 0, 0);
      glVertex3f(0, 0, 1);
      glEnd();
      }*/
    double Mx = 0;
    for (auto poly : Equestria::polygon) {
        for (auto pt : poly.pList)
            Mx = std::max(Mx, pt.len());
    }
    for (auto poly : Equestria::polygon) {
        if (poly.label >= 0) {
            auto mt = Equestria::material[poly.label].mtl;
            GLfloat Ka[] = {(float)mt.Ka.x, (float)mt.Ka.y, (float)mt.Ka.z, 1};
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, Ka);
            GLfloat Kd[] = {(float)mt.Kd.x, (float)mt.Kd.y, (float)mt.Kd.z, 1};
            glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, Kd);
            GLfloat Ks[] = {(float)mt.Ks.x, (float)mt.Ks.y, (float)mt.Ks.z, 1};
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, Ks);
            //glMaterialfv(GL_FRONT, GL_EMISSION,  sun_mat_emission);
            glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);
        }
        glBegin(GL_POLYGON);
        for (int i = 0; i < poly.num; ++i) {
            glNormal3dv(poly.normvList[i].value);
            glVertex3dv((poly.pList[i] / Mx).value);
        }
        glEnd();
    }

    glFlush();
    glutSwapBuffers();
}

void reshape(int w, int h)
{
    double rate = (double)w / h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    //glOrtho(-1, 1, -1, 1, 0, 3);
    gluPerspective(100, rate, 0.1, 3);
}

void initialize(const std::string &path)
{
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glShadeModel(GL_SMOOTH);
    //glEnable(GL_TEXTURE_2D);

    std::string spath = path;
    spath.erase(spath.rfind('/'));
    chdir(spath.c_str());

    chdir("data");
    Equestria::readModel("list.txt");
    chdir("..");

    //ptree = new Equestria::polyKDTree(Equestria::polygon.begin(), Equestria::polygon.end());
    //exit(0);
}

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
    glutInitWindowSize(WINDOW_SIZE, WINDOW_SIZE);
    glutCreateWindow("Equestrotopia");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    glutIdleFunc(display);

    initialize(argv[0]);

    glutMainLoop();
    return 0;
}
