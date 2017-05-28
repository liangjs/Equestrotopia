#include "kdtree.h"
#include "geometry.h"
#include "light.h"
#include <cmath>
#include <cstring>
#include <numeric>
#include <algorithm>

//#include <GL/glut.h>

using namespace std;

namespace Equestria {
    polyKDTree::polyKDTree(vpolyit _begin, vpolyit _end) {
        poly.assign(_begin, _end);
        vpolyit begin = poly.begin(), end = poly.end();
        for (int i = 0; i < 3; ++i)
            bdmin[i] = INF, bdmax[i] = -INF;
        for (vpolyit i = begin; i != end; ++i)
            for (int j = 0; j < 3; ++j) {
                bdmin[j] = min(bdmin[j], (*i)->bdmin[j]);
                bdmax[j] = max(bdmax[j], (*i)->bdmax[j]);
            }
        son[0] = son[1] = NULL;
        split = false;
        if (end - begin <= 1)
            return;
        vector<Polygon*> v[3];
        int besti, bestval = end - begin;
        double bestpos;
        vpolyit lend, rbegin;
        for (int i = 0; i < 3; ++i) {
            v[i].assign(begin, end);
            vpolyit _lend, _rbegin;
            double posi;
            int val = __split(v[i].begin(), v[i].end(), _lend, _rbegin, posi, i);
            if (val < bestval) {
                besti = i;
                bestval = val;
                bestpos = posi;
                lend = _lend;
                rbegin = _rbegin;
            }
        }
        if (bestval == end - begin)
            return;
        split_dir = besti, split_pos = bestpos;
        copy(v[besti].begin(), v[besti].end(), begin);
        auto rb = end - (v[besti].end() - rbegin), le = begin + (lend - v[besti].begin());
        mson = build(rb, le);
        son[0] = build(begin, rb);
        son[1] = build(le, end);
        split = true;
        //son[0] = new polyKDTree(begin, begin + (lend - v[besti].begin()));
        //son[1] = new polyKDTree(end - (v[besti].end() - rbegin), end);
    }

    polyKDTree* polyKDTree::build(vpolyit begin, vpolyit end) {
        if (begin == end)
            return NULL;
        return new polyKDTree(begin, end);
    }

    polyKDTree::~polyKDTree() {
        delete mson;
        delete son[0];
        delete son[1];
    }

    int polyKDTree::__split(vpolyit begin, vpolyit end, vpolyit& lend, vpolyit& rbegin, double& pos, int splitter) {
        int vsize = end - begin;
        vector<double> vmin(vsize), vmax(vsize), vall(vsize * 2);
        transform(begin, end, vmin.begin(), [ = ](const Polygon * p) {
            return p->bdmin[splitter];
        });
        transform(begin, end, vmax.begin(), [ = ](const Polygon * p) {
            return p->bdmax[splitter];
        });
        sort(vmin.begin(), vmin.end());
        sort(vmax.begin(), vmax.end());
        merge(vmin.begin(), vmin.end(), vmax.begin(), vmax.end(), vall.begin());
        int ans = vsize;
        pos = vall.front();
        for (auto t : vall) {
            int num_min = upper_bound(vmax.begin(), vmax.end(), t) - vmax.begin();
            int num_max = vmin.end() - lower_bound(vmin.begin(), vmin.end(), t);
            int ans_now = vsize - min(num_min, num_max);
            if (ans_now < ans) {
                ans = ans_now;
                pos = t;
            }
        }
        auto ord = [ = ](const Polygon * p) {
            return dcmp(p->bdmax[splitter], pos) <= 0 ? -1 : (dcmp(p->bdmin[splitter], pos) >= 0 ? 1 : 0);
        };
        sort(begin, end, [ & ](const Polygon * p1, const Polygon * p2) {
            return ord(p1) < ord(p2);
        });
        for (lend = begin; lend != end && ord(*lend) <= 0;)
            ++lend;
        for (rbegin = begin; rbegin != end && ord(*rbegin) < 0;)
            ++rbegin;
        return ans;
    }

    ptnKDTree::ptnKDTree(vptnit bgn, vptnit ed, int order):
        bdmin{INF, INF, INF}, bdmax{-INF, -INF, -INF}, son{NULL, NULL}, begin(bgn), end(ed) {
        for (vptnit i = bgn; i != ed; ++i)
            for (int j = 0; j < 3; ++j) {
                bdmin[j] = min(bdmin[j], (*i)->light.bgn.value[j]);
                bdmax[j] = max(bdmax[j], (*i)->light.bgn.value[j]);
            }
        int v = order;
        if (bdmin[v] == bdmax[v])
            (v += 1) %= 3;
        if (bdmin[v] == bdmax[v])
            (v += 1) %= 3;
        if (bdmin[v] == bdmax[v])
            return;
        sort(bgn, ed,
             [v](Photon * a, Photon * b)->bool{return a->light.bgn.value[v] < b->light.bgn.value[v];});
        vptnit split = bgn;
        int maxsplit = 0;
        for (vptnit i = bgn + 1; i != ed; ++i)
            if (((*i)->light.bgn.value[v] != (*(i - 1))->light.bgn.value[v]) &&
                    (maxsplit < abs((ed - i) + (i - bgn)))) {
                maxsplit = abs((ed - i) + (i - bgn));
                split = i;
            }
        son[0] = new ptnKDTree(bgn, split, (v + 1) % 3);
        son[1] = new ptnKDTree(split, ed, (v + 1) % 3);

    }

    ptnKDTree::~ptnKDTree() {
        delete son[0];
        delete son[1];
    }

    void ptnKDTree::find(const Hitpoint& hp, vector<Photon*>& lst) {
        double dx = max(abs(bdmin[0] - hp.position.x), abs(bdmax[0] - hp.position.x)),
               dy = max(abs(bdmin[1] - hp.position.y), abs(bdmax[1] - hp.position.y)),
               dz = max(abs(bdmin[2] - hp.position.z), abs(bdmax[2] - hp.position.z));
        double maxdist = sqr(dx) + sqr(dy) + sqr(dz);
        if (maxdist < hp.radius + EPS) {
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

        if (mindist > hp.radius + EPS)
            return;
        son[0]->find(hp, lst);
        son[1]->find(hp, lst);
    }
}
