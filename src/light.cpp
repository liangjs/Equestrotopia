#include "light.h"

namespace Equestria
{

    Photon::Photon(const Ray &ray, const Point &clr): light(ray), rgb(clr) {}

    Ray::Ray(): bgn(), vec() {}
    Ray::Ray(const Point &b, const Point &v): bgn(b), vec(v / v.len()) {}
}
