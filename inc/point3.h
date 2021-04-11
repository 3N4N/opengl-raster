#ifndef _POINT3_HPP_
#define _POINT3_HPP_

#include <iostream>
struct Point3 {
    double x;
    double y;
    double z;
};

std::ostream &operator<<(std::ostream &os, const Point3 &point);
std::istream &operator>>(std::istream &is, Point3 &point);

#endif //_POINT3_HPP_

