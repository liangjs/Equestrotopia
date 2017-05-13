#include "kdtree.h"
#include <cmath>
#include <cstring>
#include <numeric>
#include <algorithm>

//#include <GL/glut.h>

namespace Equestria
{
    polyKDTree::polyKDTree(vpolyit _begin, vpolyit _end): begin(_begin), end(_end)
    {
        for (int i = 0; i < 3; ++i)
            bdmin[i] = INF, bdmax[i] = -INF;
        for (vpolyit i = begin; i != end; ++i)
            for (int j = 0; j < 3; ++j) {
                bdmin[j] = std::min(bdmin[j], i->bdmin[j]);
                bdmax[j] = std::max(bdmax[j], i->bdmax[j]);
            }
        if (begin + 1 == end)
            son[0] = son[1] = NULL;
        else {
            std::vector<Polygon> v[3];
            int besti, bestval = end - begin;
            vpolyit lend, rbegin;
            for (int i = 0; i < 3; ++i) {
                v[i].assign(begin, end);
                vpolyit _lend, _rbegin;
                int val = split(v[i].begin(), v[i].end(), _lend, _rbegin, i);
                if (val < bestval) {
                    besti = i;
                    bestval = val;
                    lend = _lend;
                    rbegin = _rbegin;
                }
            }
            if (bestval == end - begin)
                son[0] = son[1] = NULL;
            else {
                std::copy(v[besti].begin(), v[besti].end(), begin);
                son[0] = new polyKDTree(begin, begin + (lend - v[besti].begin()));
                son[1] = new polyKDTree(end - (v[besti].end() - rbegin), end);
            }
        }
    }

    polyKDTree::~polyKDTree()
    {
        delete son[0];
        delete son[1];
    }

    int polyKDTree::split(vpolyit begin, vpolyit end, vpolyit &lend, vpolyit &rbegin, int splitter)
    {
        int vsize = end - begin;
        std::vector<double> vmin(vsize), vmax(vsize), vall(vsize * 2);
        std::transform(begin, end, vmin.begin(), [ = ](const Polygon & p) {return p.bdmin[splitter];});
        std::transform(begin, end, vmax.begin(), [ = ](const Polygon & p) {return p.bdmax[splitter];});
        std::sort(vmin.begin(), vmin.end());
        std::sort(vmax.begin(), vmax.end());
        std::merge(vmin.begin(), vmin.end(), vmax.begin(), vmax.end(), vall.begin());
        std::vector<double>::iterator itmin = vmin.begin(), itmax = vmax.begin();
        int ans = vsize;
        double pos = vall.front();
        for (auto t : vall) {
            int num_min = std::upper_bound(vmax.begin(), vmax.end(), t) - vmax.begin();
            int num_max = vmin.end() - std::lower_bound(vmin.begin(), vmin.end(), t);
            int ans_now = vsize - std::min(num_min, num_max);
            if (ans_now < ans) {
                ans = ans_now;
                pos = t;
            }
        }
        auto ord = [ = ](const Polygon & p) {return p.bdmax[splitter] <= pos ? -1 : (p.bdmin[splitter] >= pos ? 1 : 0);};
        std::sort(begin, end, [ & ](const Polygon & p1, const Polygon & p2) {return ord(p1) < ord(p2);});
        for (lend = begin; lend != end && ord(*lend) <= 0;)
            ++lend;
        for (rbegin = begin; rbegin != end && ord(*rbegin) < 0;)
            ++rbegin;
        return ans;
    }
    /*
        void polyKDTree::draw(double Mx)
        {
            glBegin(GL_LINE_LOOP);
            glVertex3d(bdmin[0]/Mx, bdmin[1]/Mx, bdmin[2]/Mx);
            glVertex3d(bdmin[0]/Mx, bdmin[1]/Mx, bdmax[2]/Mx);
            glVertex3d(bdmin[0]/Mx, bdmax[1]/Mx, bdmax[2]/Mx);
            glVertex3d(bdmin[0]/Mx, bdmax[1]/Mx, bdmin[2]/Mx);
            glEnd();
            glBegin(GL_LINE_LOOP);
            glVertex3d(bdmax[0]/Mx, bdmin[1]/Mx, bdmin[2]/Mx);
            glVertex3d(bdmax[0]/Mx, bdmin[1]/Mx, bdmax[2]/Mx);
            glVertex3d(bdmax[0]/Mx, bdmax[1]/Mx, bdmax[2]/Mx);
            glVertex3d(bdmax[0]/Mx, bdmax[1]/Mx, bdmin[2]/Mx);
            glEnd();
            glBegin(GL_LINES);
            glVertex3d(bdmin[0]/Mx, bdmin[1]/Mx, bdmin[2]/Mx);
            glVertex3d(bdmax[0]/Mx, bdmin[1]/Mx, bdmin[2]/Mx);
            glVertex3d(bdmin[0]/Mx, bdmin[1]/Mx, bdmax[2]/Mx);
            glVertex3d(bdmax[0]/Mx, bdmin[1]/Mx, bdmax[2]/Mx);
            glVertex3d(bdmin[0]/Mx, bdmax[1]/Mx, bdmax[2]/Mx);
            glVertex3d(bdmax[0]/Mx, bdmax[1]/Mx, bdmax[2]/Mx);
            glVertex3d(bdmin[0]/Mx, bdmax[1]/Mx, bdmin[2]/Mx);
            glVertex3d(bdmax[0]/Mx, bdmax[1]/Mx, bdmin[2]/Mx);
            glEnd();
            if (son[0])
                son[0]->draw(Mx);
            if (son[1])
                son[1]->draw(Mx);
        }
    */
    double polyKDTree::intersect(const Ray &ray, Point *p)const
    {
        double l = 0, r = INF;
        for (int i = 0; i < 3; ++i) {
            const double &di = ray.vec.value[i];
            const double &bi = ray.bgn.value[i];
            const double &mn = bdmin[i], &mx = bdmax[i];
            if (fabs(di) < EPS) {
                if (bi + EPS < mn || bi - EPS > mx)
                    return INF;
            }
            else {
                double l2 = (mn - bi) / di, r2 = (mx - bi) / di;
                if (l2 < r2)
                    l = std::max(l, l2), r = std::min(r, r2);
                else
                    l = std::max(l, r2), r = std::min(r, l2);
                if (l + EPS >= r)
                    return INF;
            }
        }
        if (son[0]) /* have 2 sons */ {
            double t = son[0]->intersect(ray, p);
            if (t != INF)
                return t;
            return son[1]->intersect(ray, p);
        }
        else {
            double ans = INF;
            for (vpolyit i = begin; i != end; ++i) {
                Point tmp;
                double t = ray.intersect(*i, &tmp);
                if (t < ans) {
                    ans = t;
                    *p = tmp;
                }
            }
            return ans;
        }
    }

    ptnKDTree::ptnKDTree(vptnit bgn, vptnit ed, int order):
        bdmin{INF, INF, INF}, bdmax{-INF, -INF, -INF}, son{NULL, NULL}, begin(bgn), end(ed)
    {
        for (vptnit i = bgn; i != ed; ++i)
            for (int j = 0; j < 3; ++j)
            {
            bdmin[j] = std::min(bdmin[j], (*i)->light.bgn.value[j]);
            bdmax[j] = std::max(bdmax[j], (*i)->light.bgn.value[j]);
            }
        int v = order;
        if (bdmin[v] == bdmax[v]) (v += 1) %= 3;
        if (bdmin[v] == bdmax[v]) (v += 1) %= 3;
        if (bdmin[v] == bdmax[v]) return;
        std::sort(bgn, ed, 
                [v](Photon *a, Photon *b)->bool{return a->light.bgn.value[v] < b->light.bgn.value[v];});
        vptnit split = bgn;
        int maxsplit = 0;
        for (vptnit i = bgn + 1; i != ed; ++i)
            if (((*i)->light.bgn.value[v] != (*(i-1))->light.bgn.value[v]) &&
                (maxsplit < abs((ed - i) + (i - bgn))))
            {
                maxsplit = abs((ed - i) + (i - bgn));
                split = i;
            }
        son[0] = new ptnKDTree(bgn, split, (v + 1) % 3);
        son[1] = new ptnKDTree(split, ed, (v + 1) % 3);

    }

    ptnKDTree::~ptnKDTree()
    {
        delete son[0];
        delete son[1];
    }

    void ptnKDTree::find(const Hitpoint &hp, std::vector<Photon*> &lst)
    {
        double dx = std::max(abs(bdmin[0] - hp.position.x), abs(bdmax[0] - hp.position.x)),
               dy = std::max(abs(bdmin[1] - hp.position.y), abs(bdmax[1] - hp.position.y)),
               dz = std::max(abs(bdmin[2] - hp.position.z), abs(bdmax[2] - hp.position.z));
        double maxdist = sqr(dx) + sqr(dy) + sqr(dz);
        if (maxdist < hp.radius + EPS)
        {
            for (vptnit i = begin; i != end; ++i)
                lst.push_back(*i);
            return;
        }

        double mindist = 0;
        for (int i = 0; i < 3; ++i)
            if (hp.position.value[i] > bdmax[i] + EPS)
                mindist += sqr(hp.position.value[i] - bdmax[i]);
            else if (hp.position.value[i] < bdmin[i] - EPS)
                mindist += sqr(hp.position.value[i] - bdmin[i]);

        if (mindist > hp.radius + EPS) return;
        son[0]->find(hp, lst);
        son[1]->find(hp, lst); 
    }
}
