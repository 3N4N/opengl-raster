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

string nm_scenefile = "test_cases/3/scene.txt";
string nm_conffile  = "test_cases/3/config.txt";
string nm_stage1file = "bin/stage1.txt";
string nm_stage2file = "bin/stage2.txt";
string nm_stage3file = "bin/stage3.txt";
string nm_zbuffile = "bin/z_buffer.txt";

int screen_width;
int screen_height;
float lim_xl;
float lim_xr;
float lim_yb;
float lim_yt;
float lim_zf;
float lim_zr;

rgb_store rainbow[6] = {
    {255,0,0},      // red
    {255,127,0},    // orange
    {255,255,0},    // yellow
    {0,255,0},      // green
    {0,0,255},      // blue
    {139,0,255}     // violet
};

typedef std::vector<std::vector<float> > vec2d;

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
    vec2d res(matA.size(), vector<float>(matB[0].size(), 0));
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
    float w;

    vec2d p_mat(4, vector<float>(1,0));
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

Point rodrigues(Point x, Point a, float angle) {

    float cost = cos(RAD(angle));
    float sint = sin(RAD(angle));


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
    float fovY, aspect_ratio, near, far;

    scenefile >> eye;
    scenefile >> look;
    scenefile >> up;
    scenefile >> fovY >> aspect_ratio >> near >> far;

    Point l = (look - eye);
    normalize(l);
    Point r = cross_product(l, up);
    normalize(r);
    Point u = cross_product(r, l);

    vec2d T(4, vector<float>(4,0));
    identity_matrix(T);
    T[0][3] = -eye.x;
    T[1][3] = -eye.y;
    T[2][3] = -eye.z;

    vec2d R(4, vector<float>(4,0));
    R[0][0] =  r.x; R[0][1] =  r.y; R[0][2] =  r.z;
    R[1][0] =  u.x; R[1][1] =  u.y; R[1][2] =  u.z;
    R[2][0] = -l.x; R[2][1] = -l.y; R[2][2] = -l.z;
    R[3][3] = 1;

    vec2d V = multiply(R, T);

    float fovX = fovY * aspect_ratio;
    float _t = near * tan(RAD(fovY / 2));
    float _r = near * tan(RAD(fovX / 2));

    vec2d P(4, vector<float>(4,0));
    P[0][0] = near/_r;
    P[1][1] = near/_t;
    P[2][2] = -(far+near)/(far-near);
    P[2][3] = -(2*far*near)/(far-near);
    P[3][2] = -1;
    // print_matrix(P);

    vec2d mat_identity(4, vector<float>(4,0));
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
            vec2d tran_mat(4, vector<float>(4,0));
            identity_matrix(tran_mat);
            scenefile >> tran_mat[0][3] >> tran_mat[1][3] >> tran_mat[2][3];
            vec2d mat = multiply(S.top(), tran_mat);
            S.pop();
            S.push(mat);
        } else if (command == "scale") {
            vec2d scale_mat(4, vector<float>(4,0));
            scenefile >> scale_mat[0][0] >> scale_mat[1][1] >> scale_mat[2][2];
            scale_mat[3][3] = 1;
            vec2d mat = multiply(S.top(), scale_mat);
            S.pop();
            S.push(mat);
        } else if (command == "rotate") {
            float angle;
            Point axis;
            scenefile >> angle >> axis;
            normalize(axis);

            Point i; i.x = 1; i.y = 0; i.z = 0;
            Point j; j.x = 0; j.y = 1; j.z = 0;
            Point k; k.x = 0; k.y = 0; k.z = 1;

            Point c1 = rodrigues(i, axis, angle);
            Point c2 = rodrigues(j, axis, angle);
            Point c3 = rodrigues(k, axis, angle);

            vec2d rot_mat(4, vector<float>(4, 0));
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
    // stagefile.open("bin/stage.txt");
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

    float dx = (lim_xr - lim_xl) / screen_width;
    float dy = (lim_yt - lim_yb) / screen_height;
    cout << dx << " " << dy << "\n";

    float top_y = lim_yt - dy/2;
    float left_x = lim_xl + dx/2;

    bitmap_image image(screen_width, screen_width);
    image.set_all_channels(0, 0, 0);

    vec2d z_buffer(screen_width, vector<float>(screen_height, lim_zr));

    float x1,x2,x3;
    float y1,y2,y3;
    float z1,z2,z3;

    for (Triangle &t : triangles) {
        x1 = t.points[0].x;
        x2 = t.points[1].x;
        x3 = t.points[2].x;
        y1 = t.points[0].y;
        y2 = t.points[1].y;
        y3 = t.points[2].y;
        z1 = t.points[0].z;
        z2 = t.points[1].z;
        z3 = t.points[2].z;

        float min_y = min(min(y1, y2), y3);
        float max_y = max(max(y1, y2), y3);
        // cout << min_y << " " << max_y <<"\n";

        int min_r = (max_y >= top_y) ? 0 : ceil((top_y-max_y)/dy);
        int max_r = (min_y <= -top_y) ? screen_height-1 : floor((top_y-min_y)/dy);
        // cout << min_r << " " << max_r <<"\n";

        float X1,X2,X3;
        float min_x,max_x;

        for (int row = min_r; row <= max_r; row++) {
            float _y = top_y - row * dy;
            if (y1-y2 == 0) {
                X1 = x2 + ((x2-x3)/(y2-y3)) * (_y - y2);
                X2 = x3 + ((x3-x1)/(y3-y1)) * (_y - y3);
            } else if (y2-y3 == 0) {
                X1 = x1 + ((x1-x2)/(y1-y2)) * (_y - y1);
                X2 = x3 + ((x3-x1)/(y3-y1)) * (_y - y3);
            } else if (y3-y1 == 0) {
                X1 = x1 + ((x1-x2)/(y1-y2)) * (_y - y1);
                X2 = x2 + ((x2-x3)/(y2-y3)) * (_y - y2);
            } else {
                X1 = x1 + ((x1-x2)/(y1-y2)) * (_y - y1);
                X2 = x2 + ((x2-x3)/(y2-y3)) * (_y - y2);
                X3 = x3 + ((x3-x1)/(y3-y1)) * (_y - y3);
                if (X1 < min(x1,min(x2,x3)) || X1 > max(x1, max(x2,x3))) {
                    swap(X1, X3);
                } else if (X2 < min(x1,min(x2,x3)) || X2 > max(x1, max(x2,x3))) {
                    swap(X2, X3);
                } else if (X3 < min(x1,min(x2,x3)) || X3 > max(x1, max(x2,x3))) {
                } else {
                    cout << "Error detected!\n";
                    exit(1);
                }
            }

            min_x = min(X1, X2);
            max_x = max(X1, X2);

            int min_c = (min_x <= left_x) ? 0 : round((min_x-left_x)/dx);
            int max_c = (max_x >= -left_x) ? screen_width-1 : round((max_x-left_x)/dx);
            // cout << min_c << " " << max_c <<"\n";

            for (int col = min_c; col <= max_c; col++) {
                float _x = left_x + col * dx;
                float a1 = x2 - x1;
                float b1 = y2 - y1;
                float c1 = z2 - z1;
                float a2 = x3 - x1;
                float b2 = y3 - y1;
                float c2 = z3 - z1;
                float a = b1 * c2 - b2 * c1;
                float b = a2 * c1 - a1 * c2;
                float c = a1 * b2 - b1 * a2;
                float d = (-a*x1 - b*y1 - c*z1);
                float Z = (-a/c)*_x + (-b/c)*_y + (-d/c);
                // cout << a << " " << b << " " << c << " " << d << " " << Z << "\n";

                if (Z > lim_zf && Z < z_buffer[row][col]) {
                    z_buffer[row][col] = Z;
                    image.set_pixel(col, row, t.color[0],t.color[1],t.color[2]);
                }
            }
            // break;
        }

        // break;
    }

    ofstream zbuffile;
    zbuffile.open(nm_zbuffile);
    if (!zbuffile.is_open()) {
        std::cerr << "Problem opening the stage3 file!\n";
        exit(1);
    }
    zbuffile << setprecision(6) << fixed;

    image.save_image("output.bmp");
    for (auto &row : z_buffer) {
        for (auto &col : row) {
            if (col < lim_zr) {
                zbuffile << col << "\t";
            }
        }
        zbuffile << "\n";
    }

    zbuffile.close();
    conffile.close();

    return 0;
}
