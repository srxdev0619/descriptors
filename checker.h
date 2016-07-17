#ifndef CHECKER_H
#define CHECKER_H
#include <common_structs.h>
#include <random>

class checker
{

    struct matrix{
        double w[3][3];
    };

    vector mul(matrix A, vector b);
    matrix Rotation(double alpha, double beta, double gamma);


public:
    checker();
    double check_rotation(int num_input_data, int num_rotations, int num_experiments, double (*func)(input_data data, double r_cut), double r_cut);
    double check_continuity(int num_input_data, int num_experiments, double epsilon_check, double r_cut, double (*func)(input_data data, double r_cut));
    input_data generate(int num, double r_cut);
    void rotate(input_data *data);
};

#endif // CHECKER_H
