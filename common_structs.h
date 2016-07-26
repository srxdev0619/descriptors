#define pi 3.14159265358979323846264338327950288419716939937510
#ifndef STRUCTS_AND_DEFINES_FOR_DESCRIPTOR3_H
#define STRUCTS_AND_DEFINES_FOR_DESCRIPTOR3_H



#include <math.h>
#include <iostream>

struct complex{
    double re, im;
    complex operator+ (complex &other);
    complex operator- (complex &other);
    complex operator* (complex &other);


    void print(){ printf("%lf %lf\n", re, im);}
};


struct vector{
    double x, y, z;
    vector operator+ (vector &other);
    vector operator- (vector &other);
    void print() {printf("%lf %lf %lf\n", x, y, z);}
};



struct input_data{
    int num_atoms;
    vector*vectors;  // from 1 to num_atoms totaly size = num_atoms + 1 including unused zero
    int *type;
    void print(FILE *f);
};

struct matrix{
    double w[3][3];
    void print();
};

vector vector_product(vector a, vector b);
vector mul(matrix A, vector b);
double distance(vector a, vector b);
matrix inverse(matrix A);
matrix createMatrix(vector a1, vector a2, vector a3);
vector create_vector(double x, double y, double z);
double abs(complex A);
double abs(vector a);
complex create_complex(double re, double im);
complex exp_if(double f);
vector mulV(vector a, double k);
double dot_product(vector a, vector b);
int max(int a, int b);
int min(int a, int b);
double minD(double a, double b);
double maxD(double a, double b);
complex conj(complex a);
double abs_double(double a);

#endif // STRUCTS_AND_DEFINES_FOR_DESCRIPTOR3_H
