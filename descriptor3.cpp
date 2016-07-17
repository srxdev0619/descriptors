#include "descriptor3.h"

descriptor3::descriptor3(int n_maximum, int l_maximum, double r_cut_in)
{
    n_max = n_maximum;
    l_max = l_maximum;
    r_cut = r_cut_in;

    q = new radial_function[n_max + 1];
    for (int i = 1; i <= n_max; i = i + 1){
        q[i].n = n_max;
        q[i].w = new double[n_max + 1];
    }


    factorial = new double[(n_max + l_max) * 10 + 11];
    factorial[0] = 1;
    for (int i = 1; i <= (n_max + l_max) * 10 + 10; i = i + 1){
        factorial[i] = factorial[i - 1] * i;
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

    lejandres = new polynomial*[l_max + 1];
    for (int i = 0; i <= l_max; i = i + 1){
        lejandres[i] = new polynomial [l_max + 1];     //lejandres[l][m] = derivative m from l polynom lejandr
    }

    for (int l = 0; l <= l_max; l = l + 1){
        for (int m = 0; m <= l_max; m = m + 1){
            lejandres[l][m] = der_m_lejandr(l, m);
        }
    }


    coef = new complex**[n_max + 1];

    for (int i = 1; i <= n_max; i = i + 1){
        coef[i] = new complex*[l_max + 1];
    }

    for (int i = 1; i <= n_max; i = i + 1){
        for (int j = 0; j <= l_max; j = j + 1){
            coef[i][j] = new complex[2 * l_max + 1];
            coef[i][j] = coef[i][j] + l_max;
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


descriptor3::polynomial descriptor3::lejandr(int n){
    polynomial ans;
    ans.n = n ;
    ans.w = new double[n + 1];
    for (int i = 0; i <= n; i = i + 1){
        ans.w[i] = 0.0;
    }

    int i0;
    if (n % 2 == 0) i0 = n / 2;
    else i0 = n/2 + 1;

    for (int i = i0; i <= n; i = i + 1){
        ans.w[2 * i - n] = factorial[2 * i] / (factorial[2 * i - n] * factorial[n-i] * factorial[i] * pow(2, n));
        if ((n-i) % 2 == 1){
            ans.w[2 * i - n] = -ans.w[2 * i - n];
        }
    }

    return ans;
}

descriptor3::polynomial descriptor3::der_m_lejandr(int n, int m){
    polynomial ans;
    ans = lejandr(n);
    for (int j = 0; j <= n - m; j = j + 1){
        ans.w[j] = ans.w[j + m] * factorial[j + m]/ factorial[j];
    }
    ans.n = n - m;
    if (ans.n < 0) {
        ans.n = 0;
        ans.w[0] = 0;
    }
    return ans;

}


double descriptor3::calculate_der_m_lejandr(int n, int m, double z){
    polynomial lej;
    lej = lejandres[n][m];
    double ans = 0;
    for (int i = 0; i <= lej.n; i = i + 1){
        ans = ans + lej.w[i] * pow(z, i);
    }
    return ans;
}




complex descriptor3::spherical_function(int l, int m, double theta, double f){
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

complex descriptor3::basis_function(int n, int l, int m, double x, double y, double z){

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

void descriptor3::process_data(input_data data){

    for (int n = 1; n <= n_max; n = n + 1){
        for (int l = 0; l <= l_max; l = l + 1){
            for (int m = -l; m <= l; m = m + 1){
                coef[n][l][m] = create_complex(0, 0);
                for (int i = 1; i <= data.num_atoms; i = i + 1){
                    if (abs(data.vectors[i]) <= r_cut){
                         complex a = basis_function(n, l, m, data.vectors[i].x, data.vectors[i].y, data.vectors[i].z);
                         coef[n][l][m] = coef[n][l][m] + a;
                    }
                }
            }
        }
    }

}

double descriptor3::calculate_power_specturm(int n, int l){
    double ans = 0;
    for (int m = -l; m <= l; m = m + 1){
        ans = ans + abs(coef[n][l][m]) * abs(coef[n][l][m]);
    }
    return ans;
}



double descriptor3::calculate_radial_function(int n, double r){

    double ans = 0;
    for (int i = 1; i <= n_max; i = i + 1){
        ans = ans + q[n].w[i] * calculate_inizial_radial_function(i, r);
    }
    return ans;
}

double descriptor3::inizial_dot_product(int alpha, int beta){
    //return pow(r_cut, alpha + beta + 6) / ((alpha + beta + 5) * (alpha + beta + 6)); // in case of *r
 //   return pow(r_cut, alpha + beta + 5) / (alpha + beta + 5); // in case of no r;
    return pow(r_cut, alpha + beta + 7) * 2 / ((alpha + beta + 5) * (alpha + beta + 6) * (alpha + beta + 7));
}


double descriptor3::dot_product(radial_function a, radial_function b){
    double sum = 0;
    for (int alpha = 1; alpha <= a.n; alpha = alpha + 1){
        for (int beta = 1; beta <= b.n; beta = beta + 1){
            sum = sum + inizial_dot_product(alpha, beta) * a.w[alpha] * b.w[beta];
        }
    }
    return sum;
}

descriptor3::radial_function descriptor3::substract(radial_function a, radial_function b){
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

descriptor3::radial_function descriptor3::mul(radial_function a, double k){
    radial_function ans;
    ans.n = a.n;
    ans.w = new double[ans.n + 1];
    for (int i = 1; i <= ans.n; i = i + 1){
        ans.w[i] = a.w[i] * k;
    }
    return ans;
}


double descriptor3::calculate_inizial_radial_function(int alpha, double r){
    return pow((r_cut - r), alpha+ 2);
}

double descriptor3::ClebshGordan(int l, int m, int l1, int m1, int l2, int m2){
    int j = l, j1 = l1, j2 = l2;
    if (m != (m1 + m2)) return 0;
    int min_possible, max_possible;
    min_possible  = abs(j1 - j2);
    max_possible = j1 + j2;
    if (j > max_possible) return 0;
    if (j < min_possible) return 0;


    double first_sqrt = sqrt((2 * j + 1) * factorial[j + j1 - j2] * factorial[j - j1 + j2] * factorial[j1 + j2 - j] / factorial[j1 + j2 + j + 1]);
    double second_sqrt = sqrt(factorial[j + m] * factorial[j - m] *
            factorial[j1 - m1] * factorial[j1 + m1]
            * factorial[j2 - m2] * factorial[j2 + m2]);

    int maxk = min(j1 + j2 - j, j1 - m1);
    maxk = min(maxk, j2 + m2);
    int mink = 0;
    mink = max(mink, j2 - j - m1);
    mink = max(mink, j1 + m2 - j);

    double sum = 0;
    double a;
    for (int k = mink; k <= maxk; k = k + 1){
        a = factorial[k] * factorial[j1 + j2 - j - k]
                * factorial[j1 - m1 - k] * factorial[j2 + m2 - k]
                * factorial[j - j2 + m1 + k] * factorial[j - j1 - m2 + k];
        a = 1.0/a;
        if (k % 2 == 1)
            a = -a;
        sum = sum + a;
    }
    return sum * first_sqrt * second_sqrt;
}

complex descriptor3::calculate_bispectrum(int n, int l, int l1, int l2){
    complex ans = create_complex(0, 0);
    complex a;
    complex k;
    for (int m = -l; m <= l; m = m + 1){
        for (int m1 = -l1; m1 <= l1; m1 = m1 + 1){
            for (int m2 = -l2; m2 <= l2; m2 = m2 + 1){\
                a = conj(coef[n][l][m]);
                a = a * coef[n][l1][m1];
                a = a * coef[n][l2][m2];
                k = create_complex(ClebshGordan(l, m, l1, m1, l2, m2), 0);
                a = a * k;
                ans = ans + a;
            }
        }
    }
    return ans;
}
