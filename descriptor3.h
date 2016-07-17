#ifndef DESCRIPTOR3_H
#define DESCRIPTOR3_H

#define pi 3.1415926
#define epsilon 0.000000001
\

#include <common_structs.h>


class descriptor3
{
private:

    struct radial_function{
        double *w;
        int n;  // that meants that it contains inizial polynomials from 1 to n inclusive (size = n + 1, including unused 0)
    };

    struct polynomial{
        double *w;
        int n; // that means that it contains coefficients from 0 to n incluseive (size = n + 1)
    };

    int n_max, l_max;
    double r_cut;
    radial_function *q;
    double *factorial;
    complex ***coef;
    polynomial **lejandres;


    polynomial lejandr(int n);
    polynomial der_m_lejandr(int n, int m);
    double calculate_der_m_lejandr(int n, int m, double z);
    complex spherical_function(int l, int m, double theta, double f);
    complex basis_function(int n, int l, int m, double x, double y, double z);
    double calculate_radial_function(int n, double r);
    double dot_product(radial_function a, radial_function b);
    radial_function substract(radial_function a, radial_function b);   
    radial_function mul(radial_function a, double k);
    double calculate_inizial_radial_function(int alpha, double r);
    double inizial_dot_product(int alpha, int beta);

public:
    double ClebshGordan(int l, int m, int l1, int m1, int l2, int m2);
    descriptor3(int n_maximum, int l_maximum, double r_cut_in);
    void process_data(input_data data);

    double calculate_power_specturm(int n, int l);
    complex calculate_bispectrum(int n, int l, int l1, int l2);
};

#endif // DESCRIPTOR3_H


