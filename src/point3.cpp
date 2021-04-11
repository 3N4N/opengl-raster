#include "point3.h"

using std::ostream;
using std::istream;

ostream &operator<<(ostream &os, const Point3 &point) {
    os << point.x << " " << point.y << " " << point.z << "\n";
    return os;
}
istream &operator>>(istream &is, Point3 &point) {
    is >> point.x >> point.y >> point.z;
    return is;
}
