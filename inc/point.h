#ifndef _POINT_H_
#define _POINT_H_

#include <iostream>
struct Point {
    double x;
    double y;
    double z;
};

std::ostream &operator<<(std::ostream &os, const Point &point);
std::istream &operator>>(std::istream &is, Point &point);
void normalize(Point &point);
double dot_product(const Point&, const Point&);
Point cross_product(const Point&, const Point&);
Point multiply(const Point&, double);
Point operator+(const Point&, const Point&);

#endif //_POINT_H_

