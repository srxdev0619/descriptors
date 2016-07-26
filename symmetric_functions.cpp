#include "symmetric_functions.h"

symmetric_functions::symmetric_functions(double r_c)
{
    r_cut = r_c;
}

double symmetric_functions::G1(input_data data){
    double ans = 0;
    double r;
    for (int i = 1; i <= data.num_atoms; i = i + 1){
        r = abs(data.vectors[i]);
        ans = ans + fc(r);
    }
    return ans;
}

double symmetric_functions::G2(input_data data, double nu, double Rs){
    double ans = 0;
    double r;
    for (int i = 1; i <= data.num_atoms; i = i + 1){
        r = abs(data.vectors[i]);
        ans = ans + fc(r) * exp(-nu * (r - Rs) * (r - Rs));
    }
    return ans;
}

double symmetric_functions::G3(input_data data, double k){
    double ans = 0;
    double r;
    for (int i = 1; i <= data.num_atoms; i = i + 1){
        r = abs(data.vectors[i]);
        ans = ans + cos(k * r) * fc(r);
    }
    return ans;
}

double symmetric_functions::G4(input_data data, double dzeta, double lambda, double nu)
{
    double ans = 0;
    vector R1, R2, R3;
    double r1, r2, r3;
    double a;
    for (int j = 1; j <= data.num_atoms; j = j + 1){
        for (int k = 1; k <= data.num_atoms; k = k + 1){
            R1 = data.vectors[j];
            R2 = data.vectors[k];
            R3 = R1 - R2;

            r1 = abs(R1);
            r2 = abs(R2);
            r3 = abs(R3);
            a = pow((1 + lambda * dot_product(R1, R2) / (r1 * r2)), dzeta)
                    * exp(-nu * (r1 * r1 + r2 * r2 + r3 * r3))
                    * fc(r1) * fc(r2) * fc(r3);
            ans = ans + a;
        }
    }
    return pow(2, 1 - dzeta) * ans;
}

double symmetric_functions::G5(input_data data, double dzeta, double lambda, double nu)
{
    double ans = 0;
    vector R1, R2, R3;
    double r1, r2;
    double a;
    for (int j = 1; j <= data.num_atoms; j = j + 1){
        for (int k = 1; k <= data.num_atoms; k = k + 1){
            R1 = data.vectors[j];
            R2 = data.vectors[k];
            R3 = R1 - R2;

            r1 = abs(R1);
            r2 = abs(R2);

            a = pow((1 + lambda * dot_product(R1, R2) / (r1 * r2)), dzeta)
                    * exp(-nu * (r1 * r1 + r2 * r2))
                    * fc(r1) * fc(r2);
            ans = ans + a;
        }
    }
    return pow(2, 1 - dzeta) * ans;
}


double symmetric_functions::fc(double r){
    if (r > r_cut) return 0;
    return 0.5 *(cos(pi * r / r_cut) + 1);
}

