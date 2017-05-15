#ifndef LIGHT_H
#define LIGHT_H

#include "geometry.h"

namespace Equestria
{
    struct Hitpoint {
        Point position;
        Point normv;
        Point raydir;
        int material;
        int x, y; // Pixel location
        Point wgt; // Pixel weight
        double radius;
        int ptncount;  // Accumulated photon count
        Point tau; // Accumulated reflected flux
        Point direct; // direct light
    };

    class Ray
    {
    public:
        Point bgn, vec;

        Ray();
        Ray(const Point &b, const Point &v);

    };

    class Photon
    {
    public:
        Ray light;
        Point rgb;

        Photon(const Ray &ray, const Point &clr);
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
        Point color;
        double power; /* 0~1 */
    };
}

#endif // LIGHT_H
