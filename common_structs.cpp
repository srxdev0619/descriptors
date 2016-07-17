#include <common_structs.h>

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
