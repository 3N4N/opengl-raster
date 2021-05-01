#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <stack>
#include <cstdlib>
#include <math.h>
#include <time.h>
#include <cstdbool>

#include "bitmap_image.h"
#include "point.h"
#include "triangle.h"

using namespace std;

// #define PI (2 * acos(0.0))
#define RAD(t) (t * (2 * acos(0.0)) / 180.0)

string nm_scenefile = "test_cases/1/scene.txt";
string nm_conffile  = "test_cases/1/config.txt";
string nm_stage1file = "bin/stage1.txt";
string nm_stage2file = "bin/stage2.txt";
string nm_stage3file = "bin/stage3.txt";

int screen_width;
int screen_height;
double lim_xl;
double lim_xr;
double lim_yb;
double lim_yt;
double lim_zf;
double lim_zr;

rgb_store rainbow[6] = {
    {255,0,0},      // red
    {255,127,0},    // orange
    {255,255,0},    // yellow
    {0,255,0},      // green
    {0,0,255},      // blue
    {139,0,255}     // violet
};

typedef std::vector<std::vector<double> > vec2d;

void print_matrix(const vec2d &mat)
{
    for (auto &row : mat) {
        for (auto &col : row) {
            cout << col << " ";
        }
        cout << "\n";
    }
}

void identity_matrix(vec2d &mat)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (i == j) mat[i][j] = 1;
            else mat[i][j] = 0;
        }
    }
}

vec2d multiply(const vec2d &matA, const vec2d &matB)
{
    vec2d res(matA.size(), vector<double>(matB[0].size(), 0));
    for (int i = 0; i < res.size(); i++) {
        for (int j = 0; j < res[0].size(); j++) {
            for (int k = 0; k < 4; k++)
                res[i][j] += matA[i][k] * matB[k][j];
        }
    }
    return res;
}

Point transformPoint(const vec2d &matrix, const Point &point)
{
    Point p;
    double w;

    vec2d p_mat(4, vector<double>(1,0));
    p_mat[0][0] = point.x;
    p_mat[1][0] = point.y;
    p_mat[2][0] = point.z;
    p_mat[3][0] = 1;

    vec2d res = multiply(matrix, p_mat);
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

    vec2d T(4, vector<double>(4,0));
    identity_matrix(T);
    T[0][3] = -eye.x;
    T[1][3] = -eye.y;
    T[2][3] = -eye.z;

    vec2d R(4, vector<double>(4,0));
    R[0][0] =  r.x; R[0][1] =  r.y; R[0][2] =  r.z;
    R[1][0] =  u.x; R[1][1] =  u.y; R[1][2] =  u.z;
    R[2][0] = -l.x; R[2][1] = -l.y; R[2][2] = -l.z;
    R[3][3] = 1;

    vec2d V = multiply(R, T);

    double fovX = fovY * aspect_ratio;
    double _t = near * tan(RAD(fovY / 2));
    double _r = near * tan(RAD(fovX / 2));

    vec2d P(4, vector<double>(4,0));
    P[0][0] = near/_r;
    P[1][1] = near/_t;
    P[2][2] = -(far+near)/(far-near);
    P[2][3] = -(2*far*near)/(far-near);
    P[3][2] = -1;
    // print_matrix(P);

    vec2d mat_identity(4, vector<double>(4,0));
    identity_matrix(mat_identity);

    string command;

    stack<vec2d> S;
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
            vec2d tran_mat(4, vector<double>(4,0));
            identity_matrix(tran_mat);
            scenefile >> tran_mat[0][3] >> tran_mat[1][3] >> tran_mat[2][3];
            vec2d mat = multiply(S.top(), tran_mat);
            S.pop();
            S.push(mat);
        } else if (command == "scale") {
            vec2d scale_mat(4, vector<double>(4,0));
            scenefile >> scale_mat[0][0] >> scale_mat[1][1] >> scale_mat[2][2];
            scale_mat[3][3] = 1;
            vec2d mat = multiply(S.top(), scale_mat);
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

            vec2d rot_mat(4, vector<double>(4, 0));
            rot_mat[0][0] = c1.x; rot_mat[0][1] = c2.x; rot_mat[0][2] = c3.x;
            rot_mat[1][0] = c1.y; rot_mat[1][1] = c2.y; rot_mat[1][2] = c3.y;
            rot_mat[2][0] = c1.z; rot_mat[2][1] = c2.z; rot_mat[2][2] = c3.z;
            rot_mat[3][3] = 1;
            vec2d mat = multiply(S.top(), rot_mat);
            S.pop();
            S.push(mat);
        } else if (command == "push") {
            vec2d new_mat = S.top();
            S.push(new_mat);
        } else if (command == "pop") {
            S.pop();
        } else if (command == "end") {
            break;
        }
    }

    scenefile.close();
    stage1file.close();
    stage2file.close();
    stage3file.close();

    ifstream conffile;
    conffile.open(nm_conffile);
    if (!conffile.is_open()) {
        std::cerr << "Problem opening the conf file!\n";
        exit(1);
    }

    conffile >> screen_width >> screen_height;
    conffile >> lim_xl >> lim_yb >> lim_zf >> lim_zr;
    lim_xr = -lim_xl;
    lim_yt = -lim_yb;


    ifstream stagefile;
    stagefile.open(nm_stage3file);
    if (!stagefile.is_open()) {
        std::cerr << "Problem opening the stage3 file!\n";
        exit(1);
    }


    Point p1, p2, p3;
    vector<Triangle> triangles;


    int c = 0;
    while(stagefile >> p1 >> p2 >> p3) {
        Triangle triangle;
        triangle.points[0] = p1;
        triangle.points[1] = p2;
        triangle.points[2] = p3;

        triangle.color[0] = rainbow[c].red;
        triangle.color[1] = rainbow[c].green;
        triangle.color[2] = rainbow[c].blue;
        c = (c + 1 == 6) ? 0 : c+1;

        triangles.push_back(triangle);
    }
    /* Provided output appears to have been made from triangles
       taken in reverse order*/
    reverse(triangles.begin(), triangles.end());

    double dx = (lim_xr - lim_xl) / screen_width;
    double dy = (lim_yt - lim_yb) / screen_height;

    double top_y = lim_yt - dy/2;
    double left_x = lim_xl + dx/2;

    bitmap_image image(screen_width, screen_width);
    image.set_all_channels(0, 0, 0);

    vec2d z_buffer(screen_width, vector<double>(screen_height, lim_zf));

    double x1,x2,x3;
    double y1,y2,y3;

    for (Triangle &t : triangles) {
        x1 = t.points[0].x;
        x2 = t.points[1].x;
        x3 = t.points[2].x;
        y1 = t.points[0].y;
        y2 = t.points[1].y;
        y3 = t.points[2].y;

        double max_y = max(max(y1, y2), y3);
        max_y = (max_y >= top_y) ? top_y : ceil(max_y/dy)*dy - dy/2;

        double min_y = min(min(y1, y2), y3);
        min_y = (min_y <= -top_y) ? -top_y : floor(min_y/dy)*dy + dy/2;

        // cout << top_y << " " << min_y << " " << max_y << " " << dy <<"\n";
        double X1,X2,X3;
        double min_x,max_x;

        for (double _y = min_y; _y < max_y; _y += dy) {
            if (y1-y2 == 0) {
                X2 = x2 + ((x2-x3)/(y2-y3)) * (_y - y2);
                X3 = x3 + ((x3-x1)/(y3-y1)) * (_y - y3);
                min_x = min(X2, X3);
                max_x = max(X2, X3);
            } else if (y2-y3 == 0) {
                X1 = x1 + ((x1-x2)/(y1-y2)) * (_y - y1);
                X3 = x3 + ((x3-x1)/(y3-y1)) * (_y - y3);
                min_x = min(X1, X3);
                max_x = max(X1, X3);
            } else if (y3-y1 == 0) {
                X1 = x1 + ((x1-x2)/(y1-y2)) * (_y - y1);
                X2 = x2 + ((x2-x3)/(y2-y3)) * (_y - y2);
                min_x = min(X1, X2);
                max_x = max(X1, X2);
            } else {
                X1 = x1 + ((x1-x2)/(y1-y2)) * (_y - y1);
                X2 = x2 + ((x2-x3)/(y2-y3)) * (_y - y2);
                X3 = x3 + ((x3-x1)/(y3-y1)) * (_y - y3);
                min_x = min(min(X1, X2), X3);
                max_x = max(max(X1, X2), X3);
            }

            min_x = (min_x <= left_x) ? left_x : floor(min_x/dx)*dx + dx/2;
            max_x = (max_x >= -left_x) ? -left_x : ceil(max_x/dx)*dx - dx/2;
            // cout << left_x << " " << min_x << " " << max_x << " " << dx <<"\n";

            for (double _x = min_x; _x < max_x; _x += dx) {
                image.set_pixel(screen_width/2+_x*screen_width/2,
                                screen_height/2-_y*screen_height/2,
                                t.color[0],t.color[1],t.color[2]);
            }
            // break;
        }

        // break;
    }

    image.save_image("output.bmp");;


    conffile.close();

    return 0;
}
