#include <QCoreApplication>
#include <descriptor3.h>
#include <checker.h>
#include <iostream>
#include <symmetric_functions.h>
#define max_check 3
int n_global, l_global, l1_global, l2_global, type_global;
double max_re, max_im;

double power_spectrum_for_check(input_data data, double r_cut){
    descriptor3 *object_descriptor3 = new descriptor3(5, 5, r_cut);
    object_descriptor3->process_data(data);
    return object_descriptor3->calculate_power_specturm(n_global, l_global);
}

double bispectrum_for_check(input_data data, double r_cut){
    descriptor3 *object_descriptor3 = new descriptor3(5, 5, r_cut);
    object_descriptor3->process_data(data);
    complex ans = object_descriptor3->calculate_bispectrum(n_global, l_global, l1_global, l2_global);
    if (abs(ans.im) > max_im)
        max_im = abs(ans.im);
    if (abs(ans.re) > max_re)
        max_re = abs(ans.re);

    if (type_global == 1) return ans.im;
    return ans.re;
}

int main(int argc, char *argv[])
{


    descriptor3 *object_descriptor3 = new descriptor3(5, 5, 1.3);
    checker *object_checker = new checker();
    input_data data = object_checker->generate(10, 1.0);

    object_descriptor3->process_data(data);
    printf("%lf", object_descriptor3->ClebshGordan(4, 2, 5, 3, 6, -1));
    FILE *f;
    f = fopen("output.txt", "w");   



    double max_error = 0;
    double a;
    for (n_global = 1; n_global <= max_check; n_global = n_global + 1){
        for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
             a = object_checker->check_rotation(10, 10, 10, power_spectrum_for_check, 1.3);
             if (a > max_error)
                 max_error = a;
        }        
    }

    fprintf(f, "max error power spectrum rotational = %lf\n", max_error);

    max_error = 0;
    for (n_global = 1; n_global <= max_check; n_global = n_global + 1){
        for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
            a = object_checker->check_continuity(10, 10, 0.001, 1.3, power_spectrum_for_check);
            if (a > max_error)
                max_error = a;
        }        
    }


    fprintf(f, "max error power spectrum continuity = %lf\n", max_error);

    max_re = 0;
    max_im = 0;

    max_error = 0;
    for (n_global = 1; n_global <= max_check; n_global = n_global + 1){
        for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
            for (l1_global = 0; l1_global <= max_check; l1_global = l1_global + 1){
                for (l2_global = 0; l2_global <= max_check; l2_global = l2_global + 1){
                    for (type_global = 0; type_global <= 1; type_global = type_global + 1){
                        a = object_checker->check_rotation(5, 3, 1, bispectrum_for_check, 1.0);
                        if (a > max_error)
                            max_error = a;
                    }
                }
            }
        }
    }



    fprintf(f, "max error bipectrum rotation = %lf\n", max_error);



   max_error = 0;
    for (n_global = 1; n_global <= max_check; n_global = n_global + 1){
        for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
            for (l1_global = 0; l1_global <= max_check; l1_global = l1_global + 1){
                for (l2_global = 0; l2_global <= max_check; l2_global = l2_global + 1){
                    for (type_global = 0; type_global <= 1; type_global = type_global + 1){
                        a = object_checker->check_continuity(10, 1, 0.001, 1.3, bispectrum_for_check);
                        if (a > max_error)
                            max_error = a;
                    }
                }
            }
        }
    }


    fprintf(f, "max error bipectrum continuity = %lf\n", max_error);
    printf("max_re = %lf\n max_im = %lf\n", max_re, max_im);

    printf("finished");
    fclose(f);
    //return a.exec();
}
