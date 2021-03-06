#include "geometry.h"
#include "brdf.h"
#include "kdtree.h"
#include "light.h"
#include "input.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <stdexcept>

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

    void Point::normalize()
    {
        double l = len();
        x /= l;
        y /= l;
        z /= l;
    }

    void Point::rotate(double dr, const Point &center, const Point &axis)
    {
        x -= center.x;
        y -= center.y;
        z -= center.z;
        double q1 = cos(dr / 2),
               q2 = sin(dr / 2) * axis.x,
               q3 = sin(dr / 2) * axis.y,
               q4 = sin(dr / 2) * axis.z;
        double tx = x, ty = y, tz = z;
        x = (sqr(q1) + sqr(q2) - sqr(q3)
             - sqr(q4)) * tx + 2 * (q2 * q3 - q1 * q4) * ty
            + 2 * (q2 * q4 + q1 * q3) * tz;
        y = 2 * (q2 * q3 + q1 * q4) * tx
            + (sqr(q1) - sqr(q2) + sqr(q3) - sqr(q4)) * ty
            + 2 * (q3 * q4 - q1 * q2) * tz;
        z = 2 * (q2 * q4 - q1 * q3) * tx
            + 2 * (q3 * q4 + q1 * q2) * ty
            + (sqr(q1) - sqr(q2) - sqr(q3) + sqr(q4)) * tz;
        x = x + center.x;
        y = y + center.y;
        z = z + center.z;
    }

    void Point::multiByChannel(const Point &p)
    {
        x *= p.x;
        y *= p.y;
        z *= p.z;
    }

    std::ostream &operator<< (std::ostream &os, const Point &p)
    {
        return os << p.x << ' ' << p.y << ' ' << p.z;
    }

    std::istream &operator>> (std::istream &is, const Point &p)
    {
        return is >> p.x >> p.y >> p.z;
    }

    Sphere::Sphere(const Point &a, double r, int label): center(a), radius(r), Object(label) {}
    //Sphere::Sphere(double ox, double oy, double oz, double r): center(ox, oy, oz), radius(r) {}

    Point Sphere::getNormal(const Point &p)const
    {
        Point v = p - center;
        return v / v.len();
    }

    void Sphere::rotate(double dr, const Point &o, const Point &axis)
    {
        center.rotate(dr, o, axis);
    }

    Polygon::Polygon(const Polygon &p): num(p.num), Object(p.label), normvf(p.normvf), pList(p.pList), normvList(p.normvList), texList(p.texList)
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
        : pList(pl), normvList(nl), texList(tl), num(pl.size()), Object(lab)
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
        Point ta = pList[c2] - pList[c1], tb = pList[c3] - pList[c1];
        xy = ta.x * tb.y - ta.y * tb.x;
        xz = ta.x * tb.z - ta.z * tb.x;
        yz = ta.y * tb.z - ta.z * tb.y;
        normvf = Point(yz, -xz, xy);
        normvf /= normvf.len();
    }

    Point Polygon::getNormal(const Point &p)const
    {
        double a0 = SignedArea(pList[c2] - p, pList[c3] - p);
        double a1 = SignedArea(pList[c3] - p, pList[c1] - p);
        double a2 = SignedArea(pList[c1] - p, pList[c2] - p);
        double s = a0 + a1 + a2;
        a0 /= s, a1 /= s, a2 /= s;
        Point N = normvList[c1] * a0 + normvList[c2] * a1 + normvList[c3] * a2;
        N /= N.len();
        return N;
    }

    void Polygon::rotate(double dr, const Point &center, const Point &axis)
    {
        for (auto &p : pList)
            p.rotate(dr, center, axis);
        for (auto &p : normvList)
            p.rotate(dr, center, axis);
        normvf.rotate(dr, center, axis);

        Point ta = pList[c2] - pList[c1], tb = pList[c3] - pList[c1];
        xy = ta.x * tb.y - ta.y * tb.x;
        xz = ta.x * tb.z - ta.z * tb.x;
        yz = ta.y * tb.z - ta.z * tb.y;

        xmin = ymin = zmin = INF;
        xmax = ymax = zmax = -INF;
        double s = 0;
        for (int i = 0; i < num; ++i) {
            xmin = std::min(xmin, pList[i].x);
            xmax = std::max(xmax, pList[i].x);
            ymin = std::min(ymin, pList[i].y);
            ymax = std::max(ymax, pList[i].y);
            zmin = std::min(zmin, pList[i].z);
            zmax = std::max(zmax, pList[i].z);
        }
    }

    void Polygon::txCoordinate(const Point &p, double &u, double &v)const
    {
        double a0 = SignedArea(pList[c2] - p, pList[c3] - p);
        double a1 = SignedArea(pList[c3] - p, pList[c1] - p);
        double a2 = SignedArea(pList[c1] - p, pList[c2] - p);
        double s = a0 + a1 + a2;
        a0 /= s, a1 /= s, a2 /= s;
        Point f = texList[c1] * a0 + texList[c2] * a1 + texList[c3] * a2;
        u = f.x, v = f.y;
    }

    double Polygon::SignedArea(const Point &p1, const Point &p2)const
    {
        return dotsProduct(crossProduct(p1, p2), normvf);
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

    Point elemMult(const Point &a, const Point &b)
    {
        return Point(a.x * b.x, a.y * b.y, a.z * b.z);
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

    double intersect(const Ray &ray, const Sphere &s, Point *p, double lasthit)
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
        if (res >= lasthit)
            return INF;
        *p = ray.bgn + res * ray.vec;
        return res;
    }

    double intersect(const Ray &ray, const Polygon &s, Point *p, double lasthit)
    {
        Point ts = ray.bgn - s.pList[s.c1];
        double k = ray.vec.x * s.yz - ray.vec.y * s.xz + ray.vec.z * s.xy,
               b = ts.x * s.yz - ts.y * s.xz + ts.z * s.xy;
        if (fabs(k) < EPS)
            return INF;
        double t = -b / k;
        if (t < EPS || t >= lasthit)
            return INF;
        Point ret = ray.bgn + t * ray.vec;

        double txy = s.normvf.x * ret.y - s.normvf.y * ret.x,
               txz = s.normvf.x * ret.z - s.normvf.z * ret.x,
               tyz = s.normvf.y * ret.z - s.normvf.z * ret.y;
        for (int i = 0; i < s.num; ++i) {
            int ni = (i + 1) % s.num;
            if (determinant(s.pList[i], s.pList[ni], s.normvf) +
                (s.pList[i].x - s.pList[ni].x) * tyz -
                (s.pList[i].y - s.pList[ni].y) * txz +
                (s.pList[i].z - s.pList[ni].z) * txy < -EPS)
                return INF;
        }

        *p = ret;
        return t;
    }

    double intersect(const Ray &ray, const polyKDTree *tree, Polygon *&p, double lasthit)
    {
        using namespace std;
        if (tree == NULL)
            return INF;
        double l = 0, r = INF;
        for (int i = 0; i < 3; ++i) {
            const double &di = ray.vec.value[i];
            const double &bi = ray.bgn.value[i];
            const double &mn = tree->bdmin[i], &mx = tree->bdmax[i];
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
                if (l - EPS >= r)
                    return INF;
            }
        }
        if (l >= lasthit)
            return INF;
        if (tree->split) {
            int k = dcmp(ray.bgn.value[tree->split_dir], tree->split_pos);
            if (k == 0) // ray.bgn.value[split_dir] == split_pos
                k = dcmp(ray.vec.value[tree->split_dir]);
            k = (k + 1) / 2; // -1 -> 0  0 -> 0   1 -> 1
            Polygon *p1, *p2, *pm;
            double t1 = intersect(ray, tree->son[k], p1, lasthit);
            double tm = intersect(ray, tree->mson, pm, min(t1, lasthit));
            bool ig = t1 != INF;
            if (tm < t1) {
                t1 = tm;
                p1 = pm;
            }
            if (ig) {
                p = p1;
                return t1;
            }
            double t2 = intersect(ray, tree->son[k ^ 1], p2, min(t1, lasthit));
            if (t1 >= lasthit && t2 >= lasthit)
                return INF;
            if (t1 > t2) {
                p = p2;
                return t2;
            }
            else {
                p = p1;
                return t1;
            }
        }
        else {
            double ans = INF;
            for (auto i = tree->poly.begin(); i != tree->poly.end(); ++i) {
                Point tmp;
                double t = intersect(ray, **i, &tmp, lasthit);
                if (t < ans) {
                    ans = t;
                    p = *i;
                }
            }
            if (ans >= lasthit)
                return INF;
            return ans;
        }
    }

    Point Material::BRDF_cos(const Point &v_in, const Point &v_out, const Point &N, double u, double v)
    {
        if (brdf == NULL) {
            Point H = v_in + v_out;
            H /= H.len();
            return mtl.getKd(u, v) * dotsProduct(v_in, N)
                + mtl.Ks * pow(std::max(0.0, dotsProduct(N, H)), mtl.Ns);
        }
        double theta_in = acos(dotsProduct(v_in, N));
        double theta_out = acos(dotsProduct(v_out, N));
        Point v_in2 = v_in - N * cos(theta_in);
        Point v_out2 = v_out - N * cos(theta_out);
        double phi = acos(dotsProduct(v_in2, v_out2) / v_in2.len() / v_out2.len());
        double r, g, b;
        BRDF::lookup_brdf_val(brdf, theta_in, phi, theta_out, 0, r, g, b);
        Point brdf = Point(r, g, b);
        if (mtl.mapKd != -1)
            brdf = (brdf + mtl.getKd(u, v)) / 2;
        return brdf * dotsProduct(v_in, N);
    }

    Ray reflect(const Point &p, const Point &N, const Point &I)
    {
        return Ray(p + N * EPS, I - 2 * dotsProduct(I, N) * N);
    }

    Ray refract(const Point &p, double n1, double n2, const Point &N, const Point &I)
    {
        double nn = n1 / n2;
        double d = dotsProduct(I, N);
        double t = 1 - nn * nn * (1 - d * d);
        if (t < 0)
            throw std::logic_error("tot inner refl");
        t = nn * d + sqrt(t);
        return Ray(p - N * EPS, nn * I - t * N);
    }

    double integral(double l, double r, const std::function<double(double)> &f)
    {
        auto simpson = [&](double a, double b) {
            double c =  a + (b - a) / 2;
            return (f(a) + 4 * f(c) + f(b)) * (b - a) / 6;
        };
        std::function<double(double, double, double)> asr = [&](double a, double b, double A) {
            double c = (b + a) / 2;
            double L = simpson(a, c), R = simpson(c, b);
            const double eps = 0.001;
            if (fabs(A - L - R) <= eps)
                return A;
            return asr(a, c, L) + asr(c, b, R);
        };
        return asr(l, r, simpson(l, r));
    }

    double intersect(const Ray &ray, Point &pos, Point &N, Object *&p)
    {
        Polygon *pl;
        double t = intersect(ray, polytree, pl);
        p = pl;
        for (auto &s : sphere) {
            Point p2;
            double t2 = intersect(ray, *s, &p2, t);
            if (t2 < t) {
                t = t2;
                p = s;
            }
        }
        if (t != INF) {
            pos = ray.bgn + t * ray.vec;
            N = p->getNormal(pos);
        }
        return t;
    }

}
