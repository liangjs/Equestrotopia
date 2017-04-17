#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>

namespace Equestria {
    const double EPS = 1E-7;
    const double INF = 1E10;

    class Point;
    class Ray;
    class Polygon;
    class Sphere;
    class AdjTable;
    class KDTree;

    class Point {
    public:
        typedef double *const double_const_ptr;
        double_const_ptr x, y, z;
        double value[3];

        Point(double x = 0, double y = 0, double z = 0);
        Point(const Point &);
        Point operator+(const Point &) const;
        Point &operator+=(const Point &);
        Point operator-(const Point &) const;
        Point &operator-=(const Point &);
        Point operator-() const;
        Point operator*(double) const;
        Point &operator*=(double);
        Point operator/(double) const;
        Point &operator/=(double);
        Point &operator=(const Point &x);

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
    private:
        double xy, xz, yz, xmin, ymin, zmin, xmax, ymax, zmax;
    public:
        int c1, c2, c3; //
        int label; // material index
        int num; // number of points
        Point normvf;
        std::vector<Point> pList, normvList, texList;

        friend class Ray;

        Polygon();
        Polygon(const std::vector<Point> &pl, const std::vector<Point> &nl,
                const std::vector<Point> &tl, int lab);
    };

    class AdjTable {
    public:
        class Node{
        public:
            int x;
            Node* nxt;

            Node();
            Node(int, Node*);
        };
        Node* head;

        AdjTable();
        void reArrange(std::vector<int>::iterator bg, std::vector<int>::iterator ed);
    };

    class KDTree {
    public:
        AdjTable lkList;
        double xmin, xmax, ymin, ymax, zmin, zmax;
        KDTree* nxt[2];

        KDTree();
        KDTree(std::vector<int>::iterator bg, std::vector<int>::iterator ed);
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
