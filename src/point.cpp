#include "point.h"

#include <cmath>

using std::ostream;
using std::istream;

ostream &operator<<(ostream &os, const Point &point)
{
    os << point.x << " " << point.y << " " << point.z << "\n";
    return os;
}

istream &operator>>(istream &is, Point &point)
{
    is >> point.x >> point.y >> point.z;
    return is;
}

void normalize(Point &point)
{
    double r = sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
    point.x /= r;
    point.y /= r;
    point.z /= r;
}

double dot_product(const Point &p1, const Point &p2)
{
    return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z;
}

Point cross_product(const Point &a, const Point &b)
{
    Point prod;

    prod.x = a.y * b.z - b.y * a.z;
    prod.y = a.z * b.x - b.z * a.x;
    prod.z = a.x * b.y - b.x * a.y;

    return prod;
}

Point multiply(const Point& p, double weight)
{
    Point q;
    q.x = p.x * weight;
    q.y = p.y * weight;
    q.z = p.z * weight;
    return q;
}

Point operator+(const Point &l, const Point &r)
{
    return {l.x + r.x, l.y + r.y, l.z + r.z};
}

Point operator-(const Point &l, const Point &r)
{
    return {l.x - r.x, l.y - r.y, l.z - r.z};
}
