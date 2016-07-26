#include "descriptor34.h"

descriptor34::descriptor34(int n_maximum, int l_maximum, double j_maximum, double r_cut_in)
{
    n_max = n_maximum;
    l_max = l_maximum;
    r_cut = r_cut_in;
    j_max = j_maximum;

    q = new radial_function[n_max + 1];
    for (int i = 1; i <= n_max; i = i + 1){
        q[i].n = n_max;
        q[i].w = new double[n_max + 1];
    }


    factorial = new double[(n_max + l_max + lrint(2 * j_maximum)) * 10 + 11];
    factorial[0] = 1;
    for (int i = 1; i <= (n_max + l_max + lrint(2 * j_maximum)) * 10 + 10; i = i + 1){
        factorial[i] = factorial[i - 1] * i;
    }


    double_factorial = new double [(n_max + l_max + lrint(2 * j_maximum)) * 10 + 11];
    double_factorial[0] = 1;
    double_factorial[1] = 1;

    for (int i = 2; i <= (n_max + l_max + lrint(2 * j_maximum) ) * 10 + 10; i = i + 1){
        double_factorial[i] = double_factorial[i - 2] * i;
    }

    for (int i = 1; i <= n_max; i = i + 1){
        for (int j = 1; j <= n_max; j = j + 1){
            if (i == j) q[i].w[j] = 1;
            else q[i].w[j] = 0;
        }
    }

    radial_function z;
    z.n = n_max;
    z.w = new double[n_max + 1];

    double k;
    for (int i = 1; i <= n_max; i = i + 1){
        for (int j = 1; j < i; j = j + 1){
            k = dot_product(q[i], q[j]) / dot_product(q[j], q[j]);
            z = mul(q[j], k);
            q[i] = substract(q[i], z);
        }
        k = dot_product(q[i], q[i]);
        q[i] = mul(q[i], 1/sqrt(k));
    }




    spherical_coef = new complex**[n_max + 1];

    for (int i = 1; i <= n_max; i = i + 1){
        spherical_coef[i] = new complex*[l_max + 1];
    }

    for (int i = 1; i <= n_max; i = i + 1){
        for (int j = 0; j <= l_max; j = j + 1){
            spherical_coef[i][j] = new complex[2 * l_max + 1];
            spherical_coef[i][j] = spherical_coef[i][j] + l_max;
        }
    }



    hyperspherical_coef = new complex**[lrint(j_maximum * 2) + 1];

    for (int i = 0; i <= lrint(j_maximum * 2); i = i + 1){
        hyperspherical_coef[i] = new complex*[2 * lrint(j_maximum * 2) + 1];
        hyperspherical_coef[i] = hyperspherical_coef[i] + lrint(j_maximum * 2);
    }


    for (int i = 0; i <= lrint(j_maximum * 2); i = i + 1){
        for (int j = -i; j <= i; j = j + 1){
            hyperspherical_coef[i][j] = new complex[2 * lrint(j_maximum * 2) + 1];
            hyperspherical_coef[i][j] = hyperspherical_coef[i][j] + lrint(j_maximum * 2);
        }
    }


    /*FILE *f = fopen("input.txt", "r");
    int num_atoms;
    fscanf(f, "%d", &num_atoms);

    vector *vectors;
    vectors = new vector[num_atoms + 1];

    for (int i = 1; i <= num_atoms; i = i + 1){
        fscanf(f, "%lf %lf %lf", &vectors[i].x, &vectors[i].y, &vectors[i].z);
    }
    fclose(f);

    FILE *g = fopen("output.txt", "w");
    calculate_and_write_power_spectrum(num_atoms, vectors, g);
    calculate_and_write_bispectrum(num_atoms, vectors, g);

    fclose(g);
    */
   /* complex co = integral(1, 1, -1,1, 1, -1);
    co.print();
    */

}

double descriptor34::smooth_function(double r){
    return (cos(pi * r / r_cut) + 1);
}

double descriptor34::calculate_gegenbauer(int alpha, int n, double z){
    double ans = 0;
    double a;
    for (int k = 0; k <= n/2; k = k + 1){
        a = factorial[n - k + alpha - 1] / (factorial[alpha - 1] * factorial[k] * factorial[n - 2 * k]);
        if (k % 2 == 1)
            a = -a;
        ans = ans + pow(2 * z, n - 2 * k) * a;
    }
    return ans;
}


complex descriptor34::hyperspherical_function(double l, double mu, double v, double theta0, double theta, double f){
    complex ans = create_complex(0, 0);
    complex a;
    complex b;
    complex k;
    for (int lambda = 0; lambda <= lrint(2 * l); lambda = lambda + 1){
        a = create_complex(0, 0);
        for (int alpha = -lambda; alpha <= lambda; alpha = alpha + 1){
            b = spherical_function(lambda, alpha, theta, f);
            k = create_complex(sqrt(4 * pi / (2 * lambda + 1)) * ClebshGordan(l, v, l, mu, lambda, alpha), 0);
            b = b * k;
            a = a + b;
        }

        k = create_complex(((double)(2 * lambda + 1)) / (2 * l + 1) * kappa(l, lambda, theta0), 0);
        a = a * k;
        switch (lambda % 4){
        case 0:
            k = create_complex(1, 0);
            break;
        case 1:
            k = create_complex(0, -1);
            break;
        case 2:
            k = create_complex(-1, 0);
            break;
        case 3:
            k = create_complex(0, 1);
            break;

        }
        a = a * k;
        ans = ans + a;
    }
    return ans;
}


double descriptor34::calculate_der_m_lejandr(int n, int m, double z){

    double ans = 0;
    int i_min;
    if ((n + m) % 2 == 0)
        i_min = (n + m)/ 2;
    else
        i_min = (n + m) / 2 + 1;
    double a;
    for (int i = i_min; i <= n; i = i + 1){
        a = factorial[2 * i] / (factorial[n - i] * factorial[i] * factorial[2 * i - n - m]);
        if ((n - i) % 2 == 1)
            a = -a;
        a = a * pow(z, 2 * i - n - m);
        ans = ans + a;
    }
    return pow(2, -n) * ans;
}




complex descriptor34::spherical_function(int l, int m, double theta, double f){
    complex a = exp_if(m * f);

    double k;
    int e;
    if (m < 0) {
        m = -m;
        e = -1;
        k = factorial[l - m] / factorial[l + m];
        if (m % 2 == 1)
            k = -k;

    }
    else{
       e = 1;
       k = 1;
    }

    double lej = calculate_der_m_lejandr(l, m, cos(theta));
    lej = lej * pow(sin(theta), m);
    lej = lej * k;
    m = m * e;

    lej = lej * sqrt((2 * l + 1) / (4 * pi) *  factorial[l - m] / factorial[l + m]);
    //if (m % 2 == 1) lej = -lej;
    complex TH = create_complex(lej, 0);
    return a * TH;
}

complex descriptor34::basis_function3(int n, int l, int m, double x, double y, double z){

    double r, theta, f;
    r = sqrt(x * x + y * y + z * z);
    if (r < epsilon){
        theta = 0;
        f=  0;
    }
    else{
        theta = acos(z / r);
        if (sqrt(x * x + y * y) < epsilon)
            f = 0;
        else{
             f = acos(x / sqrt(x * x + y * y));
             if (y < 0) f = 2 * pi - f;
        }
    }
    complex spheric = spherical_function(l, m, theta, f);
    complex b = create_complex(calculate_radial_function(n, r), 0);
    return b * spheric;
}

complex descriptor34::basis_function4(double j, double m, double m_hatch, vector position){
    double theta0, theta, f, r;
    r = abs(position);
    if (r < epsilon){
        theta = 0;
        f=  0;
    }
    else{
        theta = acos(position.z / r);
        if (sqrt(position.x * position.x + position.y * position.y) < epsilon)
            f = 0;
        else{
             f = acos(position.x / sqrt(position.x * position.x + position.y * position.y));
             if (position.y < 0) f = 2 * pi - f;
        }
    }

    theta0 = pi * abs(position) / (r_cut * r0);





    return hyperspherical_function(j, m, m_hatch, theta0, theta, f);
}

double descriptor34::kappa(double l , int lambda, double theta0){
    double a;
    a = double_factorial[2 * lambda] * sqrt(2 * l + 1) *
            sqrt(factorial[lrint(2 * l - lambda)] / factorial[lrint(2 * l + lambda + 1)])
            * pow(sin(theta0 / 2), lambda);
    return a * calculate_gegenbauer(lambda + 1, lrint(2 * l - lambda), cos(theta0/ 2));
}

void descriptor34::process_data(input_data data){
    complex k;
    for (int n = 1; n <= n_max; n = n + 1){
        for (int l = 0; l <= l_max; l = l + 1){
            for (int m = -l; m <= l; m = m + 1){
                spherical_coef[n][l][m] = create_complex(0, 0);
                for (int i = 1; i <= data.num_atoms; i = i + 1){
                    if (abs(data.vectors[i]) <= r_cut){
                         complex a = basis_function3(n, l, m, data.vectors[i].x, data.vectors[i].y, data.vectors[i].z);
                         k = create_complex(data.type[i], 0);
                         a = a * k;
                         spherical_coef[n][l][m] = spherical_coef[n][l][m] + a;
                    }
                }
            }
        }
    }

    vector zero = {0, 0, 0};
    for (int j = 0; j <= lrint(j_max * 2); j = j + 1){
        for (int m = -j; m <= j; m = m + 2){
            for (int m_hatch = -j; m_hatch <= j; m_hatch = m_hatch + 2){
                hyperspherical_coef[j][m_hatch][m] = create_complex(0, 0);
                for (int i = 1; i <= data.num_atoms; i = i + 1){
                    if (abs(data.vectors[i]) <= r_cut){
                         complex a = basis_function4(j / 2.0, m_hatch / 2.0, m / 2.0, data.vectors[i]);
                         k = create_complex(data.type[i] * smooth_function(abs(data.vectors[i])), 0);
                         a = a * k;
                         hyperspherical_coef[j][m_hatch][m] = hyperspherical_coef[j][m_hatch][m] + a;
                    }
                }

                complex a = basis_function4(j / 2.0, m_hatch / 2.0, m / 2.0, zero);
                k = create_complex(1 * smooth_function(abs(zero)), 0);
                a = a * k;
                hyperspherical_coef[j][m_hatch][m] = hyperspherical_coef[j][m_hatch][m] + a;
            }
        }
    }



}


double descriptor34::calculate_radial_function(int n, double r){

    double ans = 0;
    for (int i = 1; i <= n_max; i = i + 1){
        ans = ans + q[n].w[i] * calculate_inizial_radial_function(i, r);
    }
    return ans;
}

double descriptor34::inizial_dot_product(int alpha, int beta){
    //return pow(r_cut, alpha + beta + 6) / ((alpha + beta + 5) * (alpha + beta + 6)); // in case of *r
 //   return pow(r_cut, alpha + beta + 5) / (alpha + beta + 5); // in case of no r;
    return pow(r_cut, alpha + beta + 7) * 2 / ((alpha + beta + 5) * (alpha + beta + 6) * (alpha + beta + 7));
}


double descriptor34::dot_product(radial_function a, radial_function b){
    double sum = 0;
    for (int alpha = 1; alpha <= a.n; alpha = alpha + 1){
        for (int beta = 1; beta <= b.n; beta = beta + 1){
            sum = sum + inizial_dot_product(alpha, beta) * a.w[alpha] * b.w[beta];
        }
    }
    return sum;
}

descriptor34::radial_function descriptor34::substract(radial_function a, radial_function b){
    radial_function ans;
    ans.n = max(a.n, b.n);
    ans.w = new double[ans.n + 1];
    double a_now, b_now;
    for (int i = 1; i <= ans.n; i = i + 1){
        if (i <= a.n)
            a_now = a.w[i];
        else
            a_now = 0;

        if (i <= b.n)
            b_now = b.w[i];
        else
            b_now = 0;

        ans.w[i] = a_now - b_now;
    }
    return ans;
}

descriptor34::radial_function descriptor34::mul(radial_function a, double k){
    radial_function ans;
    ans.n = a.n;
    ans.w = new double[ans.n + 1];
    for (int i = 1; i <= ans.n; i = i + 1){
        ans.w[i] = a.w[i] * k;
    }
    return ans;
}


double descriptor34::calculate_inizial_radial_function(int alpha, double r){
    return pow((r_cut - r), alpha+ 2);
}

bool equal(double a, double b){
    if (abs_double(a - b) < epsilonClebshGordan)
        return true;
    return false;
}

bool larger(double a, double b){
    if (a - b > epsilonClebshGordan)
        return true;
    return false;
}

double descriptor34::ClebshGordan(double l, double m, double l1, double m1, double l2, double m2){

    double j = l, j1 = l1, j2 = l2;
    if (!equal(m,(m1 + m2))) return 0;
    double min_possible, max_possible;
    min_possible  = abs_double(j1 - j2);
    max_possible = j1 + j2;
    if (larger(j, max_possible)) return 0;
    if (larger(min_possible, j)) return 0;


    double first_sqrt = sqrt((2 * j + 1) * factorial[lrint(j + j1 - j2)] * factorial[lrint(j - j1 + j2)] * factorial[lrint(j1 + j2 - j)] / factorial[lrint(j1 + j2 + j + 1)]);
    double second_sqrt = sqrt(factorial[lrint(j + m)] * factorial[lrint(j - m)] *
            factorial[lrint(j1 - m1)] * factorial[lrint(j1 + m1)]
            * factorial[lrint(j2 - m2)] * factorial[lrint(j2 + m2)]);

    int maxk = min(lrint(j1 + j2 - j), lrint(j1 - m1));
    maxk = min(maxk, lrint(j2 + m2));
    int mink = 0;
    mink = max(mink, lrint(j2 - j - m1));
    mink = max(mink, lrint(j1 + m2 - j));

    double sum = 0;
    double a;
    for (int k = mink; k <= maxk; k = k + 1){
        a = factorial[k] * factorial[lrint(j1 + j2 - j) - k]
                * factorial[lrint(j1 - m1) - k] * factorial[lrint(j2 + m2) - k]
                * factorial[lrint(j - j2 + m1) + k] * factorial[lrint(j - j1 - m2) + k];
        a = 1.0/a;
        if (k % 2 == 1)
            a = -a;
        sum = sum + a;
    }
    return sum * first_sqrt * second_sqrt;
}

double descriptor34::calculate_power_specturm3(int n, int l){
    double ans = 0;
    for (int m = -l; m <= l; m = m + 1){
        ans = ans + abs(spherical_coef[n][l][m]) * abs(spherical_coef[n][l][m]);
    }
    return ans;
}

double descriptor34::calculate_bispectrum3(int n, int l, int l1, int l2){
    complex ans = create_complex(0, 0);
    complex a;
    complex k;
    for (int m = -l; m <= l; m = m + 1){
        for (int m1 = -l1; m1 <= l1; m1 = m1 + 1){
            for (int m2 = -l2; m2 <= l2; m2 = m2 + 1){\
                a = conj(spherical_coef[n][l][m]);
                a = a * spherical_coef[n][l1][m1];
                a = a * spherical_coef[n][l2][m2];
                k = create_complex(ClebshGordan(l, m, l1, m1, l2, m2), 0);
                a = a * k;
                ans = ans + a;
            }
        }
    }
    return ans.re;
}

double descriptor34::calculate_extended_power_spectrum3(int n1, int n2, int l){
    complex ans = create_complex(0, 0);
    complex a;
    for (int m = -l; m <= l; m = m + 1){
        a = conj(spherical_coef[n1][l][m]);
        a = a * spherical_coef[n2][l][m];
        ans = ans + a;
    }
    return ans.re;
}

complex descriptor34::calculate_extended_bispectrum3(int n1, int n2, int l, int l1, int l2){
    complex ans = create_complex(0, 0);
    complex a;
    complex k;
    for (int m = -l; m <= l; m = m + 1){
        for (int m1 = -l1; m1 <= l1; m1 = m1 + 1){
            for (int m2 = -l2; m2 <= l2; m2 = m2 + 1){
                a = conj(spherical_coef[n1][l][m]);
                a = a * spherical_coef[n2][l1][m1];
                a = a * spherical_coef[n2][l2][m2];
                k = create_complex(ClebshGordan(l, m, l1, m1, l2, m2), 0);
                a = a * k;
                ans = ans + a;
            }
        }
    }
    return ans;
}

 double descriptor34::calculate_bispectrum4(double j, double j1, double j2){

     complex a, b;
     complex ans;
     ans = create_complex(0.0, 0.0);

     for (int m = -lrint(2 * j); m <= lrint(2 * j); m = m + 2){
         for (int m_hatch = -lrint(2 * j); m_hatch <= lrint(2 * j); m_hatch = m_hatch + 2){
             for (int m1 = -lrint(2 * j1); m1 <= lrint(2 * j1); m1 = m1 + 2){
                 for (int m1_hatch = -lrint(2 * j1); m1_hatch <= lrint(2 * j1); m1_hatch = m1_hatch + 2){
                     for (int m2 = -lrint(2 * j2); m2 <= lrint(2 * j2); m2 = m2 + 2){
                         for (int m2_hatch = -lrint(2 * j2); m2_hatch <= lrint(2 * j2); m2_hatch = m2_hatch + 2){
                            a = create_complex(ClebshGordan(j, m / 2.0, j1, m1 / 2.0, j2, m2/ 2.0) *
                                               ClebshGordan(j, m_hatch / 2.0, j1, m1_hatch / 2.0, j2, m2_hatch / 2.0), 0);
                            a = a * hyperspherical_coef[lrint(2 * j1)][m1_hatch][m1];
                            a = a * hyperspherical_coef[lrint(2 * j2)][m2_hatch][m2];
                            b = conj(hyperspherical_coef[lrint(2 * j)][m_hatch][m]);
                            a = a * b;
                            ans = ans + a;
                         }
                     }
                 }
             }
         }
     }



     return ans.re;


 }


 double descriptor34::calculate_power_spectrum4(double j){
      double ans = 0.0;
      for (int m = -lrint(2 * j); m <= lrint(2 * j); m = m + 2){
          for (int m_hatch = -lrint(2 * j); m_hatch <= lrint(2 * j); m_hatch = m_hatch + 2){
              ans = ans + abs(hyperspherical_coef[lrint(2 * j)][m_hatch][m]) * abs(hyperspherical_coef[lrint(2 * j)][m_hatch][m]);
             // printf("%lf\n", abs(hyperspherical_coef[lrint(2 * j)][m_hatch][m]));
          }
      }
      return ans;
 }
