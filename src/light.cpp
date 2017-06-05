#include "light.h"

namespace Equestria {
    Hitpoint::Hitpoint():radius(0),ptncount(0),tau(){}

    Photon::Photon(): light(), rgb() {}
    Photon::Photon(const Ray& ray, const Point& clr, int mtl_label): light(ray), rgb(clr), mtl(mtl_label) {}
    //std::ostream& operator<< (std::ostream& os, const Photon& ptn) {
    //    return os << ptn.light << " " << ptn.rgb;
    //}

    Ray::Ray(): bgn(), vec() {}
    Ray::Ray(const Point& b, const Point& v): bgn(b), vec(v / v.len()) {}
    std::ostream& operator<< (std::ostream& os, const Ray& ray) {
        return os << ray.bgn << " " << ray.vec;
    }
}
