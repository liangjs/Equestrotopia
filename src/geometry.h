#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

const double EPS=1E-7;

namespace Equestria{

    class Point{
    public:
        double x, y, z;

        Point();
        Point(int x, int y, int z);
        Point operator+(const Point&) const;
        Point& operator+=(const Point&);
        Point operator-(const Point&) const;
        Point& operator-=(const Point&);
        Point operator-() const;
        Point operator*(double) const;
        Point& operator*=(double);
        Point operator/(double) const;
        Point& operator/=(double);

        double len() const;
        double len2() const;
    }

    class Ray{
    public:
        Point bgn, vec;

        Ray();
        Ray(const Point& b, const Point& v);

        bool intersect(const Sphere&, Point*) const;
        bool intersect(const Polygon&, Point*) const;
    }

    class Sphere{
    public:
        Point center;
        double radius;

        Sphere();
        Sphere(const Point&, double r);
        Sphere(double ox, double oy, double oz, double r);

    }

    class Polygon{
    public:
        int c1, c2, c3, label, num;
        double xy, xz, yz;
        Point normvf;
        vector<Point> pList, normvList, texList;


        Polygon();
        Polygon(const vector<Point>& pl, const vector<Point>& nl, const vector<Point>& tl, int lab, int len);
    }

    Point operator*(double, const Point&);

    double dotsProduct(const Point&, const Point&);
    Point crossProduct(const Point&, const Point&);
    double determinant(const Point&, const Point&, const Point&);
    double calcArea(const Point&, const Point&, const Point&);
    template<class T> sqr(const T& x){return x*x;};
}


#endif
