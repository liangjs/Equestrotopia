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

        typedef std::vector<Polygon>::iterator vpolyit;
        vpolyit begin, end;

        polyKDTree(vpolyit begin, vpolyit end);
        ~polyKDTree();
        //void draw(double Mx);
        double intersect(const Ray &ray, Point *p)const; /* return INF if no intersection */

        static int split(vpolyit begin, vpolyit end, vpolyit &lend, vpolyit &rbegin, int splitter);
    };

    class ptnKDTree
    {
    public:
        std::vector<Photon*> ptnlist;
        double bdmin[3], bdmax[3]; // bounding box
        ptnKDTree *son[2];

        typedef std::vector<Photon*>::iterator vptnit;
        vptnit begin, end;

        ptnKDTree(vptnit bgn, vptnit ed, int order = 0);
        ~ptnKDTree();
        void find(const Hitpoint &hp, std::vector<Photon*> &lst);
    };

    template<class T> T sqr(const T&x);
}

#endif
