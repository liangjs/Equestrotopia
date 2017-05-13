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
        double x, y; // Pixel location
        Point wgt; // Pixel weight
        double radius;
        int ptncount;  // Accumulated photon count
        Point tau; // Accumulated reflected flux
    };

    class Photon
    {
    public:
        Ray light;
        Point rgb;

        Photon(const Ray &ray, const Point &clr);
    };

}

#endif // LIGHT_H
