#ifndef OBJER_H
#define OBJER_H

/*
 * https://en.wikipedia.org/wiki/Wavefront_.obj_file
 */

#include <string>
#include <vector>
#include "geometry.h"

namespace Equestria {
    extern std::vector<Polygon> polygon; // from model.h

    void objRead(const std::string &file); // read obj file and save polygons to "polygon"
    void strSplit(const std::string &str, std::vector<std::string> &ans);
}

#endif
