#include "geometry.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

namespace Equestria {

    Point::Point(int tx, int ty, int tz): x(tx), y(ty), z(tz) {}
    Point Point::operator+(const Point &a) const {
        return Point(x + a.x, y + a.y, z + a.z);
    }
    Point &Point::operator+=(const Point &a) {
        x += a.x;
        y += a.y;
        z += a.z;
        return *this;
    }
    Point Point::operator-(const Point &a) const {
        return Point(x - a.x, y - a.y, z - a.z);
    }
    Point &Point::operator-=(const Point &a) {
        x -= a.x;
        y -= a.y;
        z -= a.z;
        return *this;
    }
    Point Point::operator-() const {
        return Point(-x, -y, -z);
    }
    Point Point::operator*(double k) const {
        return Point(k * x, k * y, k * z);
    }
    Point &Point::operator*=(double k) {
        x *= k;
        y *= k;
        z *= k;
        return *this;
    }
    Point operator*(double k, const Point &a) {
        return a * k;
    }
    Point Point::operator/(double k) const {
        return Point(x / k, y / k, z / k);
    }
    Point &Point::operator/=(double k) {
        x /= k;
        y /= k;
        z /= k;
        return *this;
    }
    double Point::len2() const {
        return sqr(x) + sqr(y) + sqr(z);
    }
    double Point::len() const {
        return sqrt(len2());
    }

    Ray::Ray(): bgn(), vec() {}
    Ray::Ray(const Point &b, const Point &v): bgn(b), vec(v / v.len()) {}

    Sphere::Sphere(): center(), radius(0) {}
    Sphere::Sphere(const Point &a, double r): center(a), radius(r) {}
    Sphere::Sphere(double ox, double oy, double oz, double r): center(ox, oy, oz), radius(r) {}

    Polygon::Polygon(): num(0), label(0) {}
    Polygon::Polygon(const std::vector<Point> &pl, const std::vector<Point> &nl,
                     const std::vector<Point> &tl, int lab)
        : pList(pl), normvList(nl), texList(tl), num(pl.size()), label(lab) {
        double s = 0;
        for (int i = 0; i < num; ++i)
            for (int j = i + 1; j < num; ++j)
                for (int k = j + 1; k < num; ++k) {
                    double ts = calcArea(pList[i], pList[j], pList[k]);
                    if (ts > s)
                        c1 = i, c2 = j, c3 = k, s = ts;
                }
        Point ta = c2 - c1, tb = c3 - c1;
        xy = ta.x * tb.y - ta.y * tb.x;
        xz = ta.x * tb.z - ta.z * tb.x;
        yz = ta.y * tb.z - ta.z * tb.y;
        normvf = Point(yz, -xz, xy);
    }

    bool Ray::intersect(const Sphere &s, Point *p) const {
        Point t = bgn - s.center;
        double b = 2 * dotsProduct(t, vec),
               c = t.len2() - sqr(s.radius),
               delta = sqr(b) - 4 * c;
        if (delta < EPS)
            return 0;
        delta = sqrt(delta);
        double t1 = (-b - delta) / 2,
               t2 = (-b + delta) / 2,
               res;
        if (t1 < EPS)
            if (t2 < EPS)
                return 0;
            else
                res = t2;
        else if (t2 < EPS)
            res = t1;
        else
            res = std::min(t1, t2);
        *p = bgn + res * vec;
        return 1;
    }

    bool Ray::intersect(const Polygon &s, Point *p) const {
        Point ts = bgn - s.pList[s.c1];
        double k = vec.x * s.yz - vec.y * s.xz + vec.z * s.xy,
               b = ts.x * s.yz - ts.y * s.xz + ts.z * s.xy;
        if (fabs(k) < EPS)
            return 0;
        double t = -b / k;
        Point ret = bgn + t * vec;

        double txy = s.normvf.x * ret.y - s.normvf.y * ret.x,
               txz = s.normvf.x * ret.z - s.normvf.z * ret.x,
               tyz = s.normvf.y * ret.z - s.normvf.z * ret.y;
        for (int i = 0; i < s.num - 1; ++i)
            if (determinant(s.pList[i], s.pList[i + 1], s.normvf) +
                    (s.pList[i].x - s.pList[i + 1].x) * tyz -
                    (s.pList[i].y - s.pList[i + 1].y) * txz +
                    (s.pList[i].z - s.pList[i + 1].z) * txy < -EPS)
                return 0;
        if (determinant(s.pList[s.num - 1], s.pList[1], s.normvf) +
                (s.pList[s.num - 1].x - s.pList[1].x) * tyz -
                (s.pList[s.num - 1].y - s.pList[1].y) * txz +
                (s.pList[s.num - 1].z - s.pList[1].z) * txy < -EPS)
            return 0;

        *p = ret;
        return 1;
    }

    double dotsProduct(const Point &a, const Point &b) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    Point crossProduct(const Point &a, const Point &b) {
        return Point(a.y * b.z - a.z * b.y,
                     -a.x * b.z + a.z * b.x,
                     a.x * b.y - a.y * b.x);
    }

    double determinant(const Point &a, const Point &b, const Point &c) {
        return a.x * b.y * c.z - a.x * b.z * c.y +
               a.y * b.z * c.x - a.y * b.x * c.z +
               a.z * b.x * c.y - a.z * b.y * c.x;
    }

    double calcArea(const Point &a, const Point &b, const Point &c) {
        return crossProduct(b - a, c - b).len() / 2;
    }


}
