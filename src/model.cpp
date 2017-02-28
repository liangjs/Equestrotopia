#include "model.h"
#include <fstream>

namespace Equestria {
    std::vector<Polygon> polygon;
    std::map<std::string, int> mtlIndex;

    void readModel(const std::string &file) {
        std::ifstream fin(file);
        fin.close();
    }
}
