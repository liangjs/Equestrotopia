#include "light.h"

namespace Equestria
{

Photon::Photon(const Ray &ray, const Point &clr): light(ray), rgb(clr) {}

}
