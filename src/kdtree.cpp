#include "kdtree.h"
#include <cstring>
#include <numeric>
#include <algorithm>

//#include <iostream>

namespace Equestria
{
    polyKDTree::polyKDTree(vpit_t _begin, vpit_t _end): begin(_begin), end(_end)
    {
        for (int i = 0; i < 3; ++i)
            bdmin[i] = INF, bdmax[i] = -INF;
        for (vpit_t i = begin; i != end; ++i)
            for (int j = 0; j < 3; ++j) {
                bdmin[j] = std::min(bdmin[j], i->bdmin[j]);
                bdmax[j] = std::max(bdmax[j], i->bdmax[j]);
            }
        if (begin + 1 == end)
            son[0] = son[1] = NULL;
        else {
            std::vector<Polygon> v[3];
            int besti, bestval = end - begin;
            vpit_t lend, rbegin;
            for (int i = 0; i < 3; ++i) {
                v[i].assign(begin, end);
                vpit_t _lend, _rbegin;
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
                //std::cout << lend - v[besti].begin() << " + " << v[besti].end() - rbegin << " ~= " << end - begin << std::endl;
                son[0] = new polyKDTree(begin, begin + (lend - v[besti].begin()));
                son[1] = new polyKDTree(end - (v[besti].end() - rbegin), end);
            }
        }
    }

    polyKDTree::~polyKDTree()
    {
        if (son[0])
            delete son[0];
        if (son[1])
            delete son[1];
    }

    int polyKDTree::split(vpit_t begin, vpit_t end, vpit_t &lend, vpit_t &rbegin, int splitter)
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

}
