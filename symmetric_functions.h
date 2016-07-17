#ifndef SYMMETRIC_FUNCTIONS_H
#define SYMMETRIC_FUNCTIONS_H

#include <common_structs.h>

class symmetric_functions
{
private:
    double r_cut;
    double fc(double r);
public:
    symmetric_functions(double r_c);
    double G1(input_data data);
    double G2(input_data data, double nu, double Rs);
    double G3(input_data data, double k);
    double G4(input_data data, double dzeta, double lambda, double nu);
    double G5(input_data data, double dzeta, double lambda, double nu);
};

#endif // SYMMETRIC_FUNCTIONS_H
