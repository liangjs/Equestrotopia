#include "geometry.h"
#include "kdtree.h"
#include "light.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>

namespace Equestria
{
    Point::Point(double tx, double ty, double tz)
    {
        value[0] = tx;
        value[1] = ty;
        value[2] = tz;
    }
    Point::Point(const Point &a)
    {
        memcpy(value, a.value, sizeof(value));
    }
    Point Point::operator+(const Point &a) const
    {
        return Point(x + a.x, y + a.y, z + a.z);
    }
    Point &Point::operator+=(const Point &a)
    {
        x += a.x;
        y += a.y;
        z += a.z;
        return *this;
    }
    Point Point::operator-(const Point &a) const
    {
        return Point(x - a.x, y - a.y, z - a.z);
    }
    Point &Point::operator-=(const Point &a)
    {
        x -= a.x;
        y -= a.y;
        z -= a.z;
        return *this;
    }
    Point Point::operator-() const
    {
        return Point(-x, -y, -z);
    }
    Point Point::operator*(double k) const
    {
        return Point(k * x, k * y, k * z);
    }
    Point &Point::operator*=(double k)
    {
        x *= k;
        y *= k;
        z *= k;
        return *this;
    }
    Point operator*(double k, const Point &a)
    {
        return a * k;
    }
    Point Point::operator/(double k) const
    {
        return Point(x / k, y / k, z / k);
    }
    Point &Point::operator/=(double k)
    {
        x /= k;
        y /= k;
        z /= k;
        return *this;
    }
    Point &Point::operator=(const Point &x)
    {
        memcpy(value, x.value, sizeof(value));
        return *this;
    }
    double Point::len2() const
    {
        return sqr(x) + sqr(y) + sqr(z);
    }
    double Point::len() const
    {
        return sqrt(len2());
    }
    std::ostream &operator<< (std::ostream &os, Point &p)
    {
        return os << p.x << ' ' << p.y << ' ' << p.z;
    }
    std::istream &operator>> (std::istream &is, Point &p)
    {
        return is >> p.x >> p.y >> p.z;
    }

    Sphere::Sphere(): center(), radius(0) {}
    Sphere::Sphere(const Point &a, double r): center(a), radius(r) {}
    Sphere::Sphere(double ox, double oy, double oz, double r): center(ox, oy, oz), radius(r) {}

    Polygon::Polygon(const Polygon &p): num(p.num), label(p.label), normvf(p.normvf), pList(p.pList), normvList(p.normvList), texList(p.texList)
    {
        memcpy(bdmin, p.bdmin, sizeof(bdmin));
        memcpy(bdmax, p.bdmax, sizeof(bdmax));
        c1 = p.c1, c2 = p.c2, c3 = p.c3;
        xy = p.xy, xz = p.xz, yz = p.yz;
    }

    Polygon &Polygon::operator= (const Polygon &p)
    {
        num = p.num, label = p.label;
        memcpy(bdmin, p.bdmin, sizeof(bdmin));
        memcpy(bdmax, p.bdmax, sizeof(bdmax));
        c1 = p.c1, c2 = p.c2, c3 = p.c3;
        normvf = p.normvf;
        xy = p.xy, xz = p.xz, yz = p.yz;
        pList = p.pList, normvList = p.normvList, texList = p.texList;
        return *this;
    }

    Polygon::Polygon(const std::vector<Point> &pl, const std::vector<Point> &nl,
                     const std::vector<Point> &tl, int lab)
        : pList(pl), normvList(nl), texList(tl), num(pl.size()), label(lab)
    {
        xmin = ymin = zmin = INF;
        xmax = ymax = zmax = -INF;
        double s = 0;
        for (int i = 0; i < num; ++i) {
            for (int j = i + 1; j < num; ++j)
                for (int k = j + 1; k < num; ++k) {
                    double ts = calcArea(pList[i], pList[j], pList[k]);
                    if (ts > s)
                        c1 = i, c2 = j, c3 = k, s = ts;
                }
            xmin = std::min(xmin, pList[i].x);
            xmax = std::max(xmax, pList[i].x);
            ymin = std::min(ymin, pList[i].y);
            ymax = std::max(ymax, pList[i].y);
            zmin = std::min(zmin, pList[i].z);
            zmax = std::max(zmax, pList[i].z);
        }
        Point ta = c2 - c1, tb = c3 - c1;
        xy = ta.x * tb.y - ta.y * tb.x;
        xz = ta.x * tb.z - ta.z * tb.x;
        yz = ta.y * tb.z - ta.z * tb.y;
        normvf = Point(yz, -xz, xy);
    }

    double dotsProduct(const Point &a, const Point &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    Point crossProduct(const Point &a, const Point &b)
    {
        return Point(a.y * b.z - a.z * b.y,
                     -a.x * b.z + a.z * b.x,
                     a.x * b.y - a.y * b.x);
    }

    double determinant(const Point &a, const Point &b, const Point &c)
    {
        return a.x * b.y * c.z - a.x * b.z * c.y +
               a.y * b.z * c.x - a.y * b.x * c.z +
               a.z * b.x * c.y - a.z * b.y * c.x;
    }

    double calcArea(const Point &a, const Point &b, const Point &c)
    {
        return crossProduct(b - a, c - b).len() / 2;
    }

    double intersect(const Ray &ray, const Sphere &s, Point *p)
    {
        Point t = ray.bgn - s.center;
        double b = 2 * dotsProduct(t, ray.vec),
               c = t.len2() - sqr(s.radius),
               delta = sqr(b) - 4 * c;
        if (delta < EPS)
            return INF;
        delta = sqrt(delta);
        double t1 = (-b - delta) / 2,
               t2 = (-b + delta) / 2,
               res;
        if (t1 < EPS)
            if (t2 < EPS)
                return INF;
            else
                res = t2;
        else if (t2 < EPS)
            res = t1;
        else
            res = std::min(t1, t2);
        *p = ray.bgn + res * ray.vec;
        return res;
    }

    double intersect(const Ray &ray, const Polygon &s, Point *p)
    {
        Point ts = ray.bgn - s.pList[s.c1];
        double k = ray.vec.x * s.yz - ray.vec.y * s.xz + ray.vec.z * s.xy,
               b = ts.x * s.yz - ts.y * s.xz + ts.z * s.xy;
        if (fabs(k) < EPS)
            return INF;
        double t = -b / k;
        if (t < EPS)
            return INF;
        Point ret = ray.bgn + t * ray.vec;

        double txy = s.normvf.x * ret.y - s.normvf.y * ret.x,
               txz = s.normvf.x * ret.z - s.normvf.z * ret.x,
               tyz = s.normvf.y * ret.z - s.normvf.z * ret.y;
        for (int i = 0; i < s.num - 1; ++i)
            if (determinant(s.pList[i], s.pList[i + 1], s.normvf) +
                (s.pList[i].x - s.pList[i + 1].x) * tyz -
                (s.pList[i].y - s.pList[i + 1].y) * txz +
                (s.pList[i].z - s.pList[i + 1].z) * txy < -EPS)
                return INF;
        if (determinant(s.pList[s.num - 1], s.pList[1], s.normvf) +
            (s.pList[s.num - 1].x - s.pList[1].x) * tyz -
            (s.pList[s.num - 1].y - s.pList[1].y) * txz +
            (s.pList[s.num - 1].z - s.pList[1].z) * txy < -EPS)
            return INF;

        *p = ret;
        return t;
    }

    double intersect(const Ray &ray, const polyKDTree *tree, Point *p)
    {
        /*double l = 0, r = INF;
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
                    l = max(l, l2), r = min(r, r2);
                else
                    l = max(l, r2), r = min(r, l2);
                if (l + EPS >= r)
                    return INF;
            }
        }*/
        typedef const double &rcd_t;
        rcd_t dx = ray.vec.x, dy = ray.vec.y, dz = ray.vec.z;
        rcd_t bx = ray.bgn.x, by = ray.bgn.y, bz = ray.bgn.z;
        rcd_t mnx = tree->bdmin[0], mny = tree->bdmin[1], mnz = tree->bdmin[2];
        rcd_t mxx = tree->bdmax[0], mxy = tree->bdmax[1], mxz = tree->bdmax[2];
        double t;
        if (bx < mnx) {
            if (dx <= 0)
                return INF;
            t = (mnx - bx) / dx;
        }
        else if (bx > mxx) {
            if (dx >= 0)
                return INF;
            t = (mxx - bx) / dx;
        }
        else if (by < mny) {
            if (dy <= 0)
                return INF;
            t = (mny - by) / dy;
        }
        else if (by > mxy) {
            if (dy >= 0)
                return INF;
            t = (mxy - by) / dy;
        }
        else if (bz < mnz) {
            if (dz <= 0)
                return INF;
            t = (mnz - bz) / dz;
        }
        else if (bz > mxz) {
            if (dz >= 0)
                return INF;
            t = (mxz - bz) / dz;
        }
        else
            goto calculate;
        {
            double tx = bx + t * dx, ty = by + t * dy, tz = bz + t * dz;
            if (tx < mnx - EPS || tx > mxx + EPS || ty < mny - EPS || ty > mxy + EPS || tz < mnz - EPS || tz > mxz + EPS)
                return INF;
        }
calculate:
        if (tree->son[0]) { /* have 2 sons */
            int k = dcmp(ray.bgn.value[tree->split_dir], tree->split_pos);
            if (k == 0) { // ray.bgn.value[split_dir] == split_pos
                k = dcmp(ray.vec.value[tree->split_dir]);
                if (k)
                    return intersect(ray, tree->son[(k + 1) / 2], p);
                else {
                    Point p0, p1;
                    double t0 = intersect(ray, tree->son[0], &p0);
                    double t1 = intersect(ray, tree->son[1], &p1);
                    if (t0 < t1) {
                        *p = p0;
                        return t0;
                    }
                    else {
                        *p = p1;
                        return t1;
                    }
                }
            }
            else {
                k = (k + 1) / 2; // -1 -> 0   1 -> 1
                double t = intersect(ray, tree->son[k], p);
                if (t != INF)
                    return t;
                return intersect(ray, tree->son[k ^ 1], p);
            }
        }
        else {
            double ans = INF;
            for (auto i = tree->begin; i != tree->end; ++i) {
                Point tmp;
                double t = intersect(ray, **i, &tmp);
                if (t < ans) {
                    ans = t;
                    *p = tmp;
                }
            }
            return ans;
        }
    }

}
