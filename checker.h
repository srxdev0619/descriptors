#ifndef CHECKER_H
#define CHECKER_H
#include "common_structs.h"


class checker
{


    matrix Rotation(double alpha, double beta, double gamma);


public:
    checker();
    double check_rotation(int num_input_data, int num_rotations, int num_experiments, double (*func)(input_data data, double r_cut), double r_cut);
    double check_continuity(int num_input_data, int num_experiments, double epsilon_check, double r_cut, double (*func)(input_data data, double r_cut));
    input_data generate(int num, double r_cut);
    void rotate(input_data *data);
    void compress(input_data *data);
    double check_compress(int num_input_data, int num_experiments, double (*func)(input_data data, double r_cut), double r_cut);
    double real_rand(double left, double right);
};

#endif // CHECKER_H
