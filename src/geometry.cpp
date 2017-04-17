#include "geometry.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>

namespace Equestria {
    Point::Point(double tx, double ty, double tz): x(value), y(value + 1), z(value + 2) {
        value[0] = tx;
        value[1] = ty;
        value[2] = tz;
    }
    Point::Point(const Point &a): x(value), y(value + 1), z(value + 2) {
        memcpy(value, a.value, sizeof(value));
    }
    Point Point::operator+(const Point &a) const {
        return Point(*x + *a.x, *y + *a.y, *z + *a.z);
    }
    Point &Point::operator+=(const Point &a) {
        *x += *a.x;
        *y += *a.y;
        *z += *a.z;
        return *this;
    }
    Point Point::operator-(const Point &a) const {
        return Point(*x - *a.x, *y - *a.y, *z - *a.z);
    }
    Point &Point::operator-=(const Point &a) {
        *x -= *a.x;
        *y -= *a.y;
        *z -= *a.z;
        return *this;
    }
    Point Point::operator-() const {
        return Point(-*x, -*y, -*z);
    }
    Point Point::operator*(double k) const {
        return Point(k * *x, k * *y, k * *z);
    }
    Point &Point::operator*=(double k) {
        *x *= k;
        *y *= k;
        *z *= k;
        return *this;
    }
    Point operator*(double k, const Point &a) {
        return a * k;
    }
    Point Point::operator/(double k) const {
        return Point(*x / k, *y / k, *z / k);
    }
    Point &Point::operator/=(double k) {
        *x /= k;
        *y /= k;
        *z /= k;
        return *this;
    }
    Point &Point::operator=(const Point &x)
    {
        memcpy(value, x.value, sizeof(value));
        return *this;
    }
    double Point::len2() const {
        return sqr(*x) + sqr(*y) + sqr(*z);
    }
    double Point::len() const {
        return sqrt(len2());
    }

    Ray::Ray(): bgn(), vec() {}
    Ray::Ray(const Point &b, const Point &v): bgn(b), vec(v / v.len()) {}

    Sphere::Sphere(): center(), radius(0) {}
    Sphere::Sphere(const Point &a, double r): center(a), radius(r) {}
    Sphere::Sphere(double ox, double oy, double oz, double r): center(ox, oy, oz), radius(r) {}

    Polygon::Polygon(): num(0), label(0), xmin(0), ymin(0), zmin(0), xmax(0), ymax(0), zmax(0) {}

    AdjTable::Node::Node(): x(0), nxt(nullptr) {}
    AdjTable::Node::Node(int loca,Node* nx): x(loca), nxt(nx) {}
    AdjTable::AdjTable(): head(nullptr) {}

    void AdjTable::reArrange(std::vector<int>::iterator bg, std::vector<int>::iterator ed) {
        Node* p = head;
        while (p) {
            Node* q = p; p = p->nxt;
            delete q;
        }
        head = nullptr;
        for (std::vector<int>::iterator i = bg; i != ed; ++i)
            head = new Node{*i, head};
    }

    KDTree::KDTree(): lkList(), nxt({nullptr, nullptr}), xmin(0), ymin(0), zmin(0), xmax(0), ymax(0), zmax(0) {}

    Polygon::Polygon(const std::vector<Point> &pl, const std::vector<Point> &nl,
                     const std::vector<Point> &tl, int lab)
        : pList(pl), normvList(nl), texList(tl), num(pl.size()), label(lab),
          xmin(INF), ymin(INF), zmin(INF), xmax(-INF), ymax(-INF), zmax(-INF) {
        double s = 0;
        for (int i = 0; i < num; ++i) {
            for (int j = i + 1; j < num; ++j)
                for (int k = j + 1; k < num; ++k) {
                    double ts = calcArea(pList[i], pList[j], pList[k]);
                    if (ts > s)
                        c1 = i, c2 = j, c3 = k, s = ts;
                }
            xmin = std::min(xmin, *pList[i].x);
            xmax = std::max(xmax, *pList[i].x);
            ymin = std::min(ymin, *pList[i].y);
            ymax = std::max(ymax, *pList[i].y);
            zmin = std::min(zmin, *pList[i].z);
            zmax = std::max(zmax, *pList[i].z);
        }
        Point ta = c2 - c1, tb = c3 - c1;
        xy = *ta.x * *tb.y - *ta.y * *tb.x;
        xz = *ta.x * *tb.z - *ta.z * *tb.x;
        yz = *ta.y * *tb.z - *ta.z * *tb.y;
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
        double k = *vec.x * s.yz - *vec.y * s.xz + *vec.z * s.xy,
               b = *ts.x * s.yz - *ts.y * s.xz + *ts.z * s.xy;
        if (fabs(k) < EPS)
            return 0;
        double t = -b / k;
        Point ret = bgn + t * vec;

        double txy = *s.normvf.x * *ret.y - *s.normvf.y * *ret.x,
               txz = *s.normvf.x * *ret.z - *s.normvf.z * *ret.x,
               tyz = *s.normvf.y * *ret.z - *s.normvf.z * *ret.y;
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
/*
    KDTree::KDTree(std::vector<int>::iterator bg, std::vector<int>::iterator ed)
            : lkList(), xmin(INF), ymin(INF), zmin(INF), xmax(-INF), ymax(-INF), zmax(-INF), nxt({nullptr, nullptr}){

        bool flag[3] = {0, 0, 0};
        int finalcnt = 0, coord = -1; double finaldivide;
        std::vector<std::pair<int, int>>  a;
        // x coordinate order
        a.clear();
        for (std::vector<std::pair<int, int>>::iterator i = bg; i != ed; ++i) {
            a.push_back(std::make_pair(polygon[*i].xmin, 1));
            a.push_back(std::make_pair(polygon[*i].xmax, -1));
            xmin = std::min(xmin, polygon[*i].xmin);
            xmax = std::max(xmax, polygon[*i].xmax);
        }
        double dividex; int cntx;
        split(a, cntx, dividex);
        if (cntx > finalcnt) coord = 0, finalcnt = cntx, finaldivide = dividex;

        // y coordinate order
        a.clear();
        for (std::vector<std::pair<int, int>>::iterator i = bg; i != ed; ++i) {
            a.push_back(std::make_pair(polygon[*i].ymin, 1));
            a.push_back(std::make_pair(polygon[*i].ymax, -1));
            ymin = std::min(ymin, polygon[*i].ymin);
            ymax = std::max(ymax, polygon[*i].ymax);
        }
        double dividey; int cnty;
        split(a, cnty, dividey);
        if (cnty > finalcnt) coord = 1, finalcnt = cnty, finaldivide = dividey;

        // z coordinate order
        a.clear();
        for (std::vector<std::pair<int, int>>::iterator i = bg; i != ed; ++i) {
            a.push_back(std::make_pair(polygon[*i].zmin, 1));
            a.push_back(std::make_pair(polygon[*i].zmax, -1));
            zmin = std::min(zmin, polygon[*i].zmin);
            zmax = std::max(zmax, polygon[*i].zmax);
        }
        double dividez; int cntz;
        split(a, cntz, dividez);
        if (cntz > finalcnt) coord = 2, finalcnt = cntz, finaldivide = dividez;

        if (!finalcnt) {
            lkList.reArrange(bg, ed);
            return;
        }

        flag[coord] = 1;
        
    }
*/
    double dotsProduct(const Point &a, const Point &b) {
        return *a.x * *b.x + *a.y * *b.y + *a.z * *b.z;
    }

    Point crossProduct(const Point &a, const Point &b) {
        return Point(*a.y * *b.z - *a.z * *b.y,
                     -*a.x * *b.z + *a.z * *b.x,
                     *a.x * *b.y - *a.y * *b.x);
    }

    double determinant(const Point &a, const Point &b, const Point &c) {
        return *a.x * *b.y * *c.z - *a.x * *b.z * *c.y +
               *a.y * *b.z * *c.x - *a.y * *b.x * *c.z +
               *a.z * *b.x * *c.y - *a.z * *b.y * *c.x;
    }

    double calcArea(const Point &a, const Point &b, const Point &c) {
        return crossProduct(b - a, c - b).len() / 2;
    }


}
