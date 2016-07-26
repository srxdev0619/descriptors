#include "common_structs.h"
const double epsilon = 0.000000000001;
const double r0 = 4.0/3.0;
const double epsilonClebshGordan = 0.001;


class descriptor34
{
private:

    struct radial_function{
        double *w;
        int n;  // that meants that it contains inizial functions from 1 to n inclusive (size = n + 1, including unused 0)
    };


    int n_max, l_max;
    double j_max;
    double r_cut;
    radial_function *q;
    double *factorial;
    double *double_factorial;
    complex ***spherical_coef;
    complex ***hyperspherical_coef;


    double calculate_der_m_lejandr(int n, int m, double z);
    double calculate_gegenbauer(int alpha, int n, double z);

    complex spherical_function(int l, int m, double theta, double f);

    complex basis_function3(int n, int l, int m, double x, double y, double z);
    complex basis_function4(double j, double m, double m_hatch, vector position);
    double calculate_radial_function(int n, double r);
    double dot_product(radial_function a, radial_function b);
    radial_function substract(radial_function a, radial_function b);   
    radial_function mul(radial_function a, double k);
    double calculate_inizial_radial_function(int alpha, double r);
    double inizial_dot_product(int alpha, int beta);
    double ClebshGordan(double l, double m, double l1, double m1, double l2, double m2);
    double kappa(double l, int lambda, double theta0);
    double smooth_function(double r);

public:
    complex hyperspherical_function(double l, double mu, double v, double theta0, double theta, double f);
    descriptor34(int n_maximum, int l_maximum, double j_maximum, double r_cut_in);
    void process_data(input_data data);
    double calculate_bispectrum4(double j, double j1, double j2);
    double calculate_power_spectrum4(double j);
    double calculate_power_specturm3(int n, int l);
    double calculate_bispectrum3(int n, int l, int l1, int l2);
    double calculate_extended_power_spectrum3(int n1, int n2, int l);
    complex calculate_extended_bispectrum3(int n1, int n2, int l, int l1, int l2);

};




