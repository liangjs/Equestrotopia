#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

const double EPS = 1E-7;

namespace Equestria {
    class Point;
    class Ray;
    class Polygon;
    class Sphere;

    class Point {
    public:
        double x, y, z;

        Point(int x = 0, int y = 0, int z = 0);
        Point operator+(const Point &) const;
        Point &operator+=(const Point &);
        Point operator-(const Point &) const;
        Point &operator-=(const Point &);
        Point operator-() const;
        Point operator*(double) const;
        Point &operator*=(double);
        Point operator/(double) const;
        Point &operator/=(double);

        double len() const;
        double len2() const;
    };

    class Ray {
    public:
        Point bgn, vec;

        Ray();
        Ray(const Point &b, const Point &v);

        bool intersect(const Sphere &, Point *) const;
        bool intersect(const Polygon &, Point *) const;
    };

    class Sphere {
    public:
        Point center;
        double radius;

        Sphere();
        Sphere(const Point &, double r);
        Sphere(double ox, double oy, double oz, double r);

    };

    class Polygon {
    public:
        int c1, c2, c3, label, num;
        double xy, xz, yz;
        Point normvf;
        std::vector<Point> pList, normvList, texList;


        Polygon();
        Polygon(const std::vector<Point> &pl, const std::vector<Point> &nl,
                const std::vector<Point> &tl, int lab);
    };

    Point operator*(double, const Point &);

    double dotsProduct(const Point &, const Point &);
    Point crossProduct(const Point &, const Point &);
    double determinant(const Point &, const Point &, const Point &);
    double calcArea(const Point &, const Point &, const Point &);
    template<class T> T sqr(const T &x) {
        return x * x;
    };
}


#endif
