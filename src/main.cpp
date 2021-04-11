#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>

#include "bitmap_image.h"
#include "point3.h"

using namespace std;

string nm_scenefile = "test_cases/1/scene.txt";

int main()
{
    ifstream scenefile;
    scenefile.open(nm_scenefile);

    if (!scenefile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }

    Point3 eye, look, up;
    double fovY, aspect_ratio, near, far;

    scenefile >> eye.x >> eye.y >> eye.z;
    scenefile >> look.x >> look.y >> look.z;
    scenefile >> up.x >> up.y >> up.z;
    scenefile >> fovY >> aspect_ratio >> near >> far;

    cout << eye.x << "," << eye.y << "," << eye.z << "\n";
    cout << look.x << "," << look.y << "," << look.z << "\n";
    cout << up.x << "," << up.y << "," << up.z << "\n";
    cout << fovY << "," << aspect_ratio << "," << near << "," << far << "\n";

    string command;

    while(true) {
        scenefile >> command;

        if (command == "triangle") {
            Point3 points[3];
            for (auto &point : points) {
                scenefile >> point;
            }
            break;
        // } else if (command == "translate") {
        // } else if (command == "scale") {
        // } else if (command == "rotate") {
        // } else if (command == "push") {
        // } else if (command == "pop") {
        // } else if (command == "end") {
        //     break;
        }
    }

    scenefile.close();

    return 0;
}
