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
string nm_stage2file = "bin/stage2.txt";
string nm_stage3file = "bin/stage3.txt";

typedef std::vector<std::vector<double> > two_d_vector;

void print_matrix(const two_d_vector &mat)
{
    for (auto &row : mat) {
        for (auto &col : row) {
            cout << col << " ";
        }
        cout << "\n";
    }
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

two_d_vector multiply(const two_d_vector &matA, const two_d_vector &matB)
{
    two_d_vector res(matA.size(), vector<double>(matB[0].size(), 0));
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[0].size(); j++) {
            for (int k = 0; k < 4; k++)
                res[i][j] += matA[i][k] * matB[k][j];
        }
    }
    return res;
}

Point transformPoint(const two_d_vector &matrix, const Point &point)
{
    Point p;
    double w;

    two_d_vector p_mat(4, vector<double>(1,0));
    p_mat[0][0] = point.x;
    p_mat[1][0] = point.y;
    p_mat[2][0] = point.z;
    p_mat[3][0] = 1;

    two_d_vector res = multiply(matrix, p_mat);
    p.x = res[0][0] / res[3][0];
    p.y = res[1][0] / res[3][0];
    p.z = res[2][0] / res[3][0];

    return p;
}

Point rodrigues(Point x, Point a, double angle) {

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

    ofstream stage2file;
    stage2file.open(nm_stage2file);
    if (!stage2file.is_open()) {
        std::cerr << "Problem opening the stage2 file!\n";
        exit(1);
    }

    ofstream stage3file;
    stage3file.open(nm_stage3file);
    if (!stage3file.is_open()) {
        std::cerr << "Problem opening the stage3 file!\n";
        exit(1);
    }

    stage1file << setprecision(7) << fixed;
    stage2file << setprecision(7) << fixed;
    stage3file << setprecision(7) << fixed;

    Point eye, look, up;
    double fovY, aspect_ratio, near, far;

    scenefile >> eye;
    scenefile >> look;
    scenefile >> up;
    scenefile >> fovY >> aspect_ratio >> near >> far;

    Point l = (look - eye);
    normalize(l);
    Point r = cross_product(l, up);
    normalize(r);
    Point u = cross_product(r, l);

    two_d_vector T(4, vector<double>(4,0));
    identity_matrix(T);
    T[0][3] = -eye.x;
    T[1][3] = -eye.y;
    T[2][3] = -eye.z;

    two_d_vector R(4, vector<double>(4,0));
    R[0][0] =  r.x; R[0][1] =  r.y; R[0][2] =  r.z;
    R[1][0] =  u.x; R[1][1] =  u.y; R[1][2] =  u.z;
    R[2][0] = -l.x; R[2][1] = -l.y; R[2][2] = -l.z;
    R[3][3] = 1;

    two_d_vector V = multiply(R, T);

    double fovX = fovY * aspect_ratio;
    double _t = near * tan(RAD(fovY / 2));
    double _r = near * tan(RAD(fovX / 2));

    two_d_vector P(4, vector<double>(4,0));
    P[0][0] = near/_r;
    P[1][1] = near/_t;
    P[2][2] = -(far+near)/(far-near);
    P[2][3] = -(2*far*near)/(far-near);
    P[3][2] = -1;
    // print_matrix(P);

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
                Point view = transformPoint(V, model);
                Point proj = transformPoint(P, view);
                stage1file << model;
                stage2file << view;
                stage3file << proj;
            }
            stage1file << "\n";
            stage2file << "\n";
            stage3file << "\n";
        } else if (command == "translate") {
            two_d_vector tran_mat(4, vector<double>(4,0));
            identity_matrix(tran_mat);
            scenefile >> tran_mat[0][3] >> tran_mat[1][3] >> tran_mat[2][3];
            two_d_vector mat = multiply(S.top(), tran_mat);
            S.pop();
            S.push(mat);
        } else if (command == "scale") {
            two_d_vector scale_mat(4, vector<double>(4,0));
            scenefile >> scale_mat[0][0] >> scale_mat[1][1] >> scale_mat[2][2];
            scale_mat[3][3] = 1;
            two_d_vector mat = multiply(S.top(), scale_mat);
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

            Point c1 = rodrigues(i, axis, angle);
            Point c2 = rodrigues(j, axis, angle);
            Point c3 = rodrigues(k, axis, angle);

            two_d_vector rot_mat(4, vector<double>(4, 0));
            rot_mat[0][0] = c1.x; rot_mat[0][1] = c2.x; rot_mat[0][2] = c3.x;
            rot_mat[1][0] = c1.y; rot_mat[1][1] = c2.y; rot_mat[1][2] = c3.y;
            rot_mat[2][0] = c1.z; rot_mat[2][1] = c2.z; rot_mat[2][2] = c3.z;
            rot_mat[3][3] = 1;
            two_d_vector mat = multiply(S.top(), rot_mat);
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
