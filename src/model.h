#ifndef MODEL_H
#define MODEL_H

#include <vector>
#include <map>
#include <string>
#include "geometry.h"

namespace Equestria {
    extern std::vector<Polygon> polygon;
    extern std::map<std::string, int> mtlIndex;

    void readModel(const std::string &file);
}

#endif
