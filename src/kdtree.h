#ifndef KDTREE_H

#include "geometry.h"
#include <vector>

namespace Equestria
{

    class polyKDTree
    {
    public:
        double bdmin[3], bdmax[3]; // boundary of 3 dimention
        polyKDTree *son[2];

        typedef std::vector<Polygon>::iterator vpit_t;
        vpit_t begin, end;

        polyKDTree(vpit_t begin, vpit_t end);
        ~polyKDTree();

        static int split(vpit_t begin, vpit_t end, vpit_t &lend, vpit_t &rbegin, int splitter);
    };

}

#endif
