#ifndef GEOMETRY_H
#define GEOMETRY_H

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
        Point bgn, vec;

        Ray();
        Ray(const Point& b, const Point& v);

        bool intersect(const Sphere&, Point*) const;
    }

    class Sphere{
    public:
        Point center;
        double radius;

        Sphere();
        Sphere(const Point&, double r);
        Sphere(double ox, double oy, double oz, double r);

    }

    Point operator*(double, const Point&);

    double dotsProduct(Point, Point);
    Point crossProduct(Point, Point);
    template<class T> sqr(const T& x){return x*x;};
}


#endif
