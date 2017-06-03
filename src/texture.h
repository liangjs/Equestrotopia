#ifndef TEXTURE_H
#define TEXTURE_H

#include "geometry.h"

namespace Equestria
{
    class Texture
    {
    private:
        Point get_pixel(int x, int y);
    public:
        int w, h;
        unsigned char *buf;
        Point get(double u, double v);
    };
}
#endif
