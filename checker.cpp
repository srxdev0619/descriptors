#include "checker.h"
#include <time.h>

checker::checker()
{
    srand(time(0));
}

double real_rand(double left, double right){
    return ((right - left)* rand()) / RAND_MAX + left;
}

input_data checker::generate(int num, double r_cut){
    input_data ans;
    ans.num_atoms = num;
    ans.vectors = new vector[ans.num_atoms + 1];

    for (int i = 1; i <= num; i = i + 1){
          ans.vectors[i].x = real_rand(-r_cut, r_cut);
          ans.vectors[i].y = real_rand(-r_cut, r_cut);
          ans.vectors[i].z = real_rand(-r_cut, r_cut);

    }
    return ans;
}

checker::matrix checker::Rotation(double alpha, double beta, double gamma){
        matrix A;
        A.w[0][0] = cos(alpha) * cos(gamma) - sin(alpha) * cos(beta) * sin(gamma);
        A.w[1][0] = sin(alpha) * cos(gamma) + cos(alpha) * cos(beta) * sin(gamma);
        A.w[2][0] = sin(beta) * sin(gamma);

        A.w[0][1] = -cos(alpha) * sin(gamma) - sin(alpha) * cos(beta) * cos(gamma);
        A.w[1][1] = -sin(alpha) * sin(gamma) + cos(alpha) * cos(beta) * cos(gamma);
        A.w[2][1] = sin(beta) * cos(gamma);

        A.w[0][2] = sin(alpha) * sin(beta);
        A.w[1][2] = -cos(alpha) * sin(beta);
        A.w[2][2] = cos(beta);
        return A;


}


vector checker::mul(matrix A, vector b){
    vector ans;
    ans.x = A.w[0][0] * b.x + A.w[0][1] * b.y + A.w[0][2] * b.z;
    ans.y = A.w[1][0] * b.x + A.w[1][1] * b.y + A.w[1][2] * b.z;
    ans.z = A.w[2][0] * b.x + A.w[2][1] * b.y + A.w[2][2] * b.z;
    return ans;
}

void checker::rotate(input_data *data){
    double alpha = real_rand(0, 10);
    double beta = real_rand(0, 10);
    double gamma = real_rand(0, 10);
    matrix R = Rotation(alpha, beta, gamma);
    for (int i = 1; i <= data->num_atoms; i = i + 1){
        data->vectors[i] = mul(R, data->vectors[i]);
    }
}


double checker::check_rotation(int num_input_data, int num_rotations, int num_experiments, double (*func)(input_data data, double r_cut), double r_cut){

    double max_error = 0;
    for (int i = 1; i <= num_experiments; i = i + 1){
        input_data data = generate(num_input_data, r_cut);
        double min = func(data, r_cut);
        double max = func(data, r_cut);
        double now;
        for (int j = 1; j <= num_rotations; j = j + 1){
           rotate(&data);
           now = func(data, r_cut);
           if (now < min) min = now;
           if (now > max) max = now;
        }

        if (max - min > max_error)
            max_error = max - min;

    }

    return max_error;

}

double abs(double a){
    if (a < 0) return -a;
    return a;
}

double checker::check_continuity(int num_input_data, int num_experiments, double epsilon_check, double r_cut, double (*func)(input_data data, double r_cut)){
    double max_error = 0;
    for (int i = 1; i <= num_experiments; i = i + 1){
        input_data data = generate(num_input_data, r_cut);

        for (int j = 1; j <= num_input_data; j = j + 1){
            data.vectors[j] = mulV(data.vectors[j], (1.0 - epsilon_check) / abs(data.vectors[j]) * r_cut);
            double a = func(data, r_cut);
            data.vectors[j] = mulV(data.vectors[j], (1.0 + epsilon_check) / abs(data.vectors[j]) * r_cut);
            double b = func(data, r_cut);
            if (abs(a - b) > max_error)
                max_error = abs(a - b);
        }
    }
    return max_error;
}
