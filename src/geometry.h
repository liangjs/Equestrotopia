#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <vector>
#include <iostream>
#include <functional>

namespace Equestria {
    const double EPS = 1E-7;
    const double INF = 1E10;

    class Point;
    class Polygon;
    class Sphere;
    class Ray;
    class polyKDTree;

    typedef double& double_ref;

    class Point {
        public:
            double value[3];
            double_ref x = value[0],
                       y = value[1],
                       z = value[2];

            Point(double x = 0, double y = 0, double z = 0);
            Point(const Point&);
            Point operator+(const Point&) const;
            Point& operator+=(const Point&);
            Point operator-(const Point&) const;
            Point& operator-=(const Point&);
            Point operator-() const;
            Point operator*(double) const;
            Point& operator*=(double);
            Point operator/(double) const;
            Point& operator/=(double);
            Point& operator=(const Point& x);

            double len() const;
            double len2() const;
            void normalize();
            void rotate(double dr, const Point& axis);
            void multiByChannel(const Point& p);
            // http://www.zhihu.com/question/23005815
            friend std::ostream& operator<< (std::ostream& os, const Point& p);
            friend std::istream& operator>> (std::istream& is, const Point& p);
    };

    class Sphere {
        public:
            Point center;
            double radius;

            Sphere();
            Sphere(const Point&, double r);
            Sphere(double ox, double oy, double oz, double r);

    };

    class Polygon {
        public:
            double xy, xz, yz; // used to optimize calculation

            double bdmin[3], bdmax[3];
            double_ref xmin = bdmin[0], ymin = bdmin[1], zmin = bdmin[2], xmax = bdmax[0], ymax = bdmax[1], zmax = bdmax[2];
            int c1, c2, c3; // c1,c2,c3 consist of triangle in polygon with maximum area
            int label; // material index
            int num; // number of points
            Point normvf;
            std::vector<Point> pList, normvList, texList;

            Polygon(const Polygon& p);
            Polygon(const std::vector<Point>& pl, const std::vector<Point>& nl,
                    const std::vector<Point>& tl, int lab);
            Polygon& operator= (const Polygon& p);
            Point getNormal(const Point& p)const;
            void rotate(double dr, const Point& axis);
            void txCoordinate(const Point &p, double &u, double &v)const;
            double SignedArea(const Point &p1, const Point &p2)const;
    };

    Point operator*(double, const Point&);

    double dotsProduct(const Point&, const Point&);
    Point crossProduct(const Point&, const Point&);
    Point elemMult(const Point&, const Point&);
    double determinant(const Point&, const Point&, const Point&);
    double calcArea(const Point&, const Point&, const Point&);
    template<class T> inline T sqr(const T& x) {
        return x * x;
    };
    inline int dcmp(double x, double y = 0) {
        return fabs(x - y) < EPS ? 0 : (x < y ? -1 : 1);
    }

    /* return INF if no intersection */
    double intersect(const Ray &ray, const Sphere &, Point *);
    double intersect(const Ray &ray, const Polygon &, Point *,double lasthit=INF);
    double intersect(const Ray &ray, const polyKDTree *, Polygon *&p,double lasthit=INF);
    
    Ray reflect(const Point &p, const Point &N, const Point &I);
    Ray refract(const Point &p, double n1, double n2, const Point &N, const Point &I);
}


#endif
