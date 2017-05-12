#ifndef KDTREE_H

#include "geometry.h"
#include <vector>

namespace Equestria
{

    class polyKDTree
    {
    public:
        double bdmin[3], bdmax[3]; // bounding box
        polyKDTree *son[2];

        typedef std::vector<Polygon>::iterator vpit_t;
        vpit_t begin, end;

        polyKDTree(vpit_t begin, vpit_t end);
        ~polyKDTree();
        //void draw(double Mx);
        double intersect(const Ray &ray, Point *p)const; /* return INF if no intersection */

        static int split(vpit_t begin, vpit_t end, vpit_t &lend, vpit_t &rbegin, int splitter);
    };

    class ptnKDTree
    {
    public:
        std::vector<Photon*> ptnlist;
        double bdmin[3], bdmax[3]; // bounding box
        ptnKDTree *son[2];

        ptnKDTree(int order, std::vector<Photon*>::iterator bgn, std::vector<Photon*>::iterator ed);
        ~ptnKDTree();
    };

}

#endif
