#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stack>
#include <cstdlib>

#include "bitmap_image.h"
#include "point.h"

using namespace std;

// #define PI (2 * acos(0.0))
#define RAD(t) (t * (2 * acos(0.0)) / 180.0)

string nm_scenefile = "test_cases/1/scene.txt";
string nm_stage1file = "bin/stage1.txt";

typedef std::vector<std::vector<double> > two_d_vector;

void print_matrix(two_d_vector &mat)
{
    cout << mat[0][0] <<" "<< mat[0][1] <<" "<< mat[0][2] <<" "<< mat[0][3] << "\n";
    cout << mat[1][0] <<" "<< mat[1][1] <<" "<< mat[1][2] <<" "<< mat[1][3] << "\n";
    cout << mat[2][0] <<" "<< mat[2][1] <<" "<< mat[2][2] <<" "<< mat[2][3] << "\n";
    cout << mat[3][0] <<" "<< mat[3][1] <<" "<< mat[3][2] <<" "<< mat[3][3] << "\n";
}

void identity_matrix(two_d_vector &mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) mat[i][j] = 1;
            else mat[i][j] = 0;
        }
    }
}

Point transformPoint(const two_d_vector &matrix, const Point &point)
{
    Point p;
    double w;

    w = matrix[3][0]*point.x + matrix[3][1]*point.y + matrix[3][2]*point.z + matrix[3][3];
    p.x = matrix[0][0]*point.x + matrix[0][1]*point.y + matrix[0][2]*point.z + matrix[0][3];
    p.y = matrix[1][0]*point.x + matrix[1][1]*point.y + matrix[1][2]*point.z + matrix[1][3];
    p.z = matrix[2][0]*point.x + matrix[2][1]*point.y + matrix[2][2]*point.z + matrix[2][3];

    return p;
}

two_d_vector multiply_matrix( const two_d_vector &mat1,
                              const two_d_vector &mat2)
{
    two_d_vector res(4, vector<double>(4, 0));
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++)
                res[i][j] += mat1[i][k] * mat2[k][j];
        }
    }
    return res;
}

Point formulaOfRodrigues(Point x, Point a, double angle) {

    double cost = cos(RAD(angle));
    double sint = sin(RAD(angle));


    return multiply(x, cost)
        + multiply(a, (1-cost)*dot_product(a,x))
        + multiply(cross_product(a,x), sint);
}

int main()
{
    ifstream scenefile;
    scenefile.open(nm_scenefile);

    if (!scenefile.is_open()) {
        std::cerr << "Problem opening the scene file!\n";
        exit(1);
    }

    ofstream stage1file;
    stage1file.open(nm_stage1file);

    if (!stage1file.is_open()) {
        std::cerr << "Problem opening the stage1 file!\n";
        exit(1);
    }

    stage1file << setprecision(3) << fixed;

    Point eye, look, up;
    double fovY, aspect_ratio, near, far;

    scenefile >> eye;
    scenefile >> look;
    scenefile >> up;
    scenefile >> fovY >> aspect_ratio >> near >> far;

    two_d_vector mat_identity(4, vector<double>(4,0));
    identity_matrix(mat_identity);

    string command;

    stack<two_d_vector> S;
    S.push(mat_identity);

    while(true) {
        scenefile >> command;

        if (command == "triangle") {
            Point points[3];
            for (auto &point : points) {
                scenefile >> point;
                Point model = transformPoint(S.top(), point);
                stage1file << model;
            }
            stage1file << "\n";
        } else if (command == "translate") {
            two_d_vector tran_mat(4, vector<double>(4,0));
            identity_matrix(tran_mat);
            scenefile >> tran_mat[0][3] >> tran_mat[1][3] >> tran_mat[2][3];
            two_d_vector mat = multiply_matrix(S.top(), tran_mat);
            S.pop();
            S.push(mat);
        } else if (command == "scale") {
            two_d_vector scale_mat(4, vector<double>(4,0));
            scenefile >> scale_mat[0][0] >> scale_mat[1][1] >> scale_mat[2][2];
            scale_mat[3][3] = 1;
            two_d_vector mat = multiply_matrix(S.top(), scale_mat);
            S.pop();
            S.push(mat);
        } else if (command == "rotate") {
            double angle;
            Point axis;
            scenefile >> angle >> axis;
            normalize(axis);

            Point i; i.x = 1; i.y = 0; i.z = 0;
            Point j; j.x = 0; j.y = 1; j.z = 0;
            Point k; k.x = 0; k.y = 0; k.z = 1;

            Point c1 = formulaOfRodrigues(i, axis, angle);
            Point c2 = formulaOfRodrigues(j, axis, angle);
            Point c3 = formulaOfRodrigues(k, axis, angle);

            two_d_vector rot_mat(4, vector<double>(4, 0));
            rot_mat[0][0] = c1.x; rot_mat[0][1] = c2.x; rot_mat[0][2] = c3.x;
            rot_mat[1][0] = c1.y; rot_mat[1][1] = c2.y; rot_mat[1][2] = c3.y;
            rot_mat[2][0] = c1.z; rot_mat[2][1] = c2.z; rot_mat[2][2] = c3.z;
            rot_mat[3][3] = 1;
            two_d_vector mat = multiply_matrix(S.top(), rot_mat);
            S.pop();
            S.push(mat);
        } else if (command == "push") {
            two_d_vector new_mat = S.top();
            S.push(new_mat);
        } else if (command == "pop") {
            S.pop();
        } else if (command == "end") {
            break;
        }
    }

    scenefile.close();

    return 0;
}
