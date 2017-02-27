#include "geometry.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

namespace Equestria{

    Point::Point():x(0),y(0),z(0){}
    Point::Point(int tx, int ty, int tz):x(tx),y(ty),z(tz){}
    Point Point::operator+(const Point& a) const{return Point(x + a.x, y + a.y, z + a.z);}
    Point& Point::operator+=(const Point &a){x += a.x; y += a.y; z += a.z; return *this;}
    Point Point::operator-(const Point &a) const{return Point(x - a.x, y - a.y, z - a.z);}
    Point& Point::operator-=(const Point &a){x -= a.x; y -= a.y; z -= a.z; return *this;}
    Point Point::operator-() const{return Point(-x, -y, -z);}
    Point Point::operator*(double k) const{return Point(k * x, k * y, k * z);}
    Point& Point::operator*=(double k){x *= k; y *= k; z *= k; return *this;}
    Point Point::operator/(double k) const{return Point(x / k, y / k, z / k);}
    Point& Point::operator/=(double k){x /= k; y /= k; z /= k; return *this;}
    double Point::len2() const{return sqr(x)+sqr(y)+sqr(z);}
    double Point::len() const{return sqrt(len2());}

    Ray::Ray():bgn(),vec(){}
    Ray::Ray(const Point& b, const Point& v):bgn(b),vec(v/v.len()){}

    Sphere::Sphere():center(),radius(0){}
    Sphere::Sphere(const Point& a, double r):center(a),radius(r){}
    Sphere::Sphere(double ox, double oy, double oz, double r):center(ox, oy, oz),radius(r){}

    bool Ray::intersect(const Sphere& s, Point* p) const{
        Point t = bgn - s.center;
        double b = 2 * dotsProduct(t, vec),
               c = t.len2() - sqr(s.radius),
               delta = sqr(b) - 4 * c;
        if (delta < EPS) return 0;
        delta = sqrt(delta);
        double t1 = (-b - delta) / 2,
               t2 = (-b + delta) / 2,
               res;
        if (t1 < EPS) if (t2 < EPS) return 0;
                         else res = t2;
           else if (t2 < EPS) res = t1;
                   else res = std::min(t1, t2);
        *p = bgn + res * vec;
        return 1;
    }



    Point operator*(double k, const Point& a){return a*k}

    double dotsProduct(Point a, Point b){return a.x * b.x + a.y * b.y + a.z * b.z;}
    Point crossProduct(Point a, Point b){
        return Point(a.y * b.z - a.z * b.y,
                    -a.x * b.z + a.z * b.x,
                     a.x * b.y - a.y * b.x);
    }


}
