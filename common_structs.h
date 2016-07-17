#ifndef STRUCTS_AND_DEFINES_FOR_DESCRIPTOR3_H
#define STRUCTS_AND_DEFINES_FOR_DESCRIPTOR3_H

#define pi 3.1415926

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
};



struct input_data{
    int num_atoms;
    vector*vectors;  // from 1 to num_atoms totaly size = num_atoms + 1 including unused zero
    void print(FILE *f);
};

vector create_vector(double x, double y, double z);
double abs(complex A);
double abs(vector a);
complex create_complex(double re, double im);
complex exp_if(double f);
vector mulV(vector a, double k);
double dot_product(vector a, vector b);
int max(int a, int b);
int min(int a, int b);
complex conj(complex a);

#endif // STRUCTS_AND_DEFINES_FOR_DESCRIPTOR3_H
