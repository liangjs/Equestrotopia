#include "texture.h"
#include <cmath>

namespace Equestria
{
    Point Texture::get_pixel(int x, int y)
    {
        unsigned char *ptr = buf + 3 * (x * w + y);
        unsigned char r, g, b;
        r = *ptr, g = *(ptr + 1), b = *(ptr + 2);
        return Point(r, g, b) / 255;
    }

    Point Texture::get(double u, double v)
    {
        if (u < 0)
            u = 0;
        if (v < 0)
            v = 0;
        if (u > 1)
            u = 1;
        if (v > 1)
            v = 1;
        double x = (1 - v) * (h - 1), y = u * (w - 1);
        int x0 = floor(x), y0 = floor(y);
        if (x0 == h - 1 || y0 == w - 1)
            return get_pixel(x0, y0);
        int x1 = x0 + 1, y1 = y0 + 1;
        double a = x - x0, b = y - y0;
        return (1 - a) * (1 - b) * get_pixel(x0, y0) + a * (1 - b) * get_pixel(x1, y0) + (1 - a) * b * get_pixel(x0, y1) + a * b * get_pixel(x1, y1);
    }
}
