#ifndef LIGHT_H
#define LIGHT_H

#include <iostream>
#include "geometry.h"

namespace Equestria {
    struct Hitpoint {
        Point position;
        Point normv;
        Point raydir;
        int mtl_label;
        double u, v;
        int x, y; // Pixel location
        Point wgt; // Pixel weight
        double radius;
        double ptncount;  // Accumulated photon count
        Point tau; // Accumulated reflected flux
        Point direct; // direct light

        Hitpoint();
    };

    class Ray {
        public:
            Point bgn, vec;

            Ray();
            Ray(const Point& b, const Point& v);
            friend std::ostream& operator<< (std::ostream&, const Ray&);
            inline void Print(FILE *file) { bgn.Print(file); vec.Print(file); }
    };

    class Photon {
        public:
            Ray light;
            Point rgb;
            int mtl;

            Photon();
            Photon(const Ray& ray, const Point& clr, int mtl_label);
            //friend std::ostream& operator<< (std::ostream&, const Photon&);
    };

    struct Camera {
        Point focus;
        Point o, vx, vy;
        Point normal; // normal = (o - focus) / |o - focus|
        // |o - focus| = focal length
        // o ~ (height/2,width/2)
        // o-vx/2->o+vx/2 ~ (0,width/2)->(height,width/2)
        // o-vy/2->o+vy/2 ~ (height/2,0)->(height/2,width)
        /*
                |
                |
         -------o--------> y
                |
                |
                |
                |
                v
                x
         */
    };

    struct Light {
        Point pos;
        Point I; /* Radiant Intensity */
        double power;
    };
}

#endif // LIGHT_H
