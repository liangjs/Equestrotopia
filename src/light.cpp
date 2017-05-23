#include "light.h"

namespace Equestria {

    Photon::Photon(const Ray& ray, const Point& clr): light(ray), rgb(clr) {}
    std::ostream& operator<< (std::ostream& os, const Photon& ptn) {
        return os << ptn.light << " color: " << ptn.rgb;
    }

    Ray::Ray(): bgn(), vec() {}
    Ray::Ray(const Point& b, const Point& v): bgn(b), vec(v / v.len()) {}
    std::ostream& operator<< (std::ostream& os, const Ray& ray) {
        return os << "ray_bgn: " << ray.bgn << " ray_vec: " << ray.vec;
    }
}
