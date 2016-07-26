#include "common_structs.h"

complex complex::operator+(complex &other){
    return create_complex(re + other.re, im + other.im);
}

complex conj(complex a){
    return create_complex(a.re, -a.im);
}

complex complex::operator-(complex &other){
    return create_complex(re - other.re, im - other.im);
}

complex complex::operator*(complex &other){
    return create_complex(re * other.re - im * other.im, re * other.im + im * other.re);
}

vector create_vector(double x, double y, double z){
    vector ans;
    ans.x = x;
    ans.y = y;
    ans.z = z;
    return ans;
}

vector vector::operator+ (vector &other){
    return create_vector(x + other.x, y + other.y, z + other.z);
}

vector vector::operator- (vector &other){
    return create_vector(x - other.x, y - other.y, z - other.z);
}

double abs(complex A){
    return sqrt(A.re * A.re + A.im * A.im);
}

complex create_complex(double re, double im){
    complex ans;
    ans.re =re;
    ans.im = im;
    return ans;
}

complex exp_if(double f){
    return create_complex(cos(f), sin(f));
}

vector mulV(vector a, double k){
    vector ans;
    ans.x = a.x * k;
    ans.y = a.y * k;
    ans.z = a.z * k;
    return ans;
}

vector mul(matrix A, vector b){
    vector ans;
    ans.x = A.w[0][0] * b.x + A.w[0][1] * b.y + A.w[0][2] * b.z;
    ans.y = A.w[1][0] * b.x + A.w[1][1] * b.y + A.w[1][2] * b.z;
    ans.z = A.w[2][0] * b.x + A.w[2][1] * b.y + A.w[2][2] * b.z;
    return ans;
}



double abs(vector a){
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

void input_data::print(FILE *f){
    fprintf(f, "%d\n", num_atoms);
    for (int i = 1; i <= num_atoms; i = i + 1){
        fprintf(f, "%lf %lf %lf\n", vectors[i].x, vectors[i].y, vectors[i].z);
    }

}

double dot_product(vector a, vector b){
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}


int max(int a, int b){
    if (a > b) return a;
    return b;
}

int min(int a, int b){
    if (a < b) return a;
    return b;
}

double minD(double a, double b){
    if (a < b) return a;
    return b;
}

double maxD(double a, double b){
    if (a > b) return a;
    return b;
}

double abs_double(double a){
    if (a < 0) return -a;
    return a;
}
matrix createMatrix(vector a1, vector a2, vector a3){
    matrix A;

    A.w[0][0] = a1.x;
    A.w[1][0] = a1.y;
    A.w[2][0] = a1.z;

    A.w[0][1] = a2.x;
    A.w[1][1] = a2.y;
    A.w[2][1] = a2.z;

    A.w[0][2] = a3.x;
    A.w[1][2] = a3.y;
    A.w[2][2] = a3.z;
    return A;

}

matrix inverse(matrix A){
    double determinant =    A.w[0][0]*(A.w[1][1]*A.w[2][2]-A.w[2][1]*A.w[1][2])
                            -A.w[0][1]*(A.w[1][0]*A.w[2][2]-A.w[1][2]*A.w[2][0])
                            +A.w[0][2]*(A.w[1][0]*A.w[2][1]-A.w[1][1]*A.w[2][0]);
    double invdet = 1/determinant;
    matrix ans;
    ans.w[0][0] =  (A.w[1][1]*A.w[2][2]-A.w[2][1]*A.w[1][2])*invdet;
    ans.w[0][1] = -(A.w[0][1]*A.w[2][2]-A.w[0][2]*A.w[2][1])*invdet;
    ans.w[0][2] =  (A.w[0][1]*A.w[1][2]-A.w[0][2]*A.w[1][1])*invdet;
    ans.w[1][0] = -(A.w[1][0]*A.w[2][2]-A.w[1][2]*A.w[2][0])*invdet;
    ans.w[1][1] =  (A.w[0][0]*A.w[2][2]-A.w[0][2]*A.w[2][0])*invdet;
    ans.w[1][2] = -(A.w[0][0]*A.w[1][2]-A.w[1][0]*A.w[0][2])*invdet;
    ans.w[2][0] =  (A.w[1][0]*A.w[2][1]-A.w[2][0]*A.w[1][1])*invdet;
    ans.w[2][1] = -(A.w[0][0]*A.w[2][1]-A.w[2][0]*A.w[0][1])*invdet;
    ans.w[2][2] =  (A.w[0][0]*A.w[1][1]-A.w[1][0]*A.w[0][1])*invdet;
    return ans;
}


vector vector_product(vector a, vector b){
    vector ans;
    ans.x = a.y * b.z - a.z * b.y;
    ans.y = a.z * b.x - a.x * b.z;
    ans.z = a.x * b.y - a.y * b.x;
    return ans;
}

void matrix::print(){
    for (int i = 0; i < 3; i = i + 1){
        for (int j = 0; j < 3; j = j + 1){
            printf("%lf ", w[i][j]);
        }
        printf("\n");
    }
}
double distance(vector a, vector b){
    return abs(a - b);
}
