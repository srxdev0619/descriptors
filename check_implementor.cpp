#include "check_implementor.h"
#include "checker.h"

int n_global, l_global, l1_global, l2_global, n1_global,  n2_global;
int type_global;
double max_global;
double max_re_global, max_im_global;

double power_spectrum_for_check(input_data data, double r_cut){
    descriptor34 *object_descriptor34 = new descriptor34(max_check, max_check, 0,  r_cut);
    object_descriptor34->process_data(data);
    double ans = object_descriptor34->calculate_power_specturm3(n_global, l_global);
    if (abs_double(ans) > max_global)
        max_global = abs_double(ans);
    return ans;
}

double bispectrum_for_check(input_data data, double r_cut){
    descriptor34 *object_descriptor34 = new descriptor34(max_check, max_check, 0,  r_cut);
    object_descriptor34->process_data(data);
    double ans = object_descriptor34->calculate_bispectrum3(n_global, l_global, l1_global, l2_global);
    if (abs_double(ans) > max_global)
        max_global = abs_double(ans);

    return ans;
}

double bispectrum4_for_check(input_data data, double r_cut){
    descriptor34 *object_descriptor34 = new descriptor34(0, 0, max_check / 2.0, r_cut);
    object_descriptor34->process_data(data);
    double ans = object_descriptor34->calculate_bispectrum4(n_global / 2.0, n1_global / 2.0, n2_global/ 2.0);
    if (abs_double(ans) > max_global)
        max_global = abs_double(ans);

    return ans;
}

double power_spectrum4_for_check(input_data data, double r_cut){
    descriptor34 *object_descriptor34 = new descriptor34(0, 0, max_check / 2.0, r_cut);
    object_descriptor34->process_data(data);
    double ans = object_descriptor34->calculate_power_spectrum4(n_global / 2.0);
    if (abs_double(ans) > max_global)
        max_global = abs_double(ans);
    return ans;
}

double extended_power_spectrum_for_check(input_data data, double r_cut){
    descriptor34 *object_descriptor34 = new descriptor34(max_check, max_check, 0,  r_cut);
    object_descriptor34->process_data(data);
    double ans = object_descriptor34->calculate_extended_power_spectrum3(n1_global, n2_global, l_global);
    if (abs_double(ans) > max_global)
        max_global = abs_double(ans);
    return ans;
}
double extended_bispectrum_for_check(input_data data, double r_cut){
    descriptor34 *object_descriptor34 = new descriptor34(max_check, max_check, 0, r_cut);
    object_descriptor34->process_data(data);
    complex ans = object_descriptor34->calculate_extended_bispectrum3(n1_global, n2_global, l_global, l1_global, l2_global);

    if (abs_double(ans.re) > max_re_global)
       max_re_global = abs_double(ans.re);
    if (abs_double(ans.im) > max_im_global)
        max_im_global = abs_double(ans.im);
    if (type_global == 0)
         return ans.re;
    return ans.im;
}



void check_implementor::check_all(){

    FILE *f;
    f = fopen("output.txt", "w");


    max_global = 0;
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
            a = object_checker->check_continuity(10, 10, epsilon_continuity, 1.3, power_spectrum_for_check);
            if (a > max_error)
                max_error = a;
        }
    }


    fprintf(f, "max error power spectrum continuity = %lf\n", max_error);
    fprintf(f, "max abs of power spectrum = %lf\n\n", max_global);

    max_global = 0;
    max_error = 0;
    for (n_global = 1; n_global <= max_check; n_global = n_global + 1){
        for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
            for (l1_global = 0; l1_global <= max_check; l1_global = l1_global + 1){
                for (l2_global = 0; l2_global <= max_check; l2_global = l2_global + 1){
                    a = object_checker->check_rotation(5, 3, 1, bispectrum_for_check, 1.3);
                    if (a > max_error)
                        max_error = a;

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
                   a = object_checker->check_continuity(10, 1, epsilon_continuity, 1.3, bispectrum_for_check);
                   if (a > max_error)
                      max_error = a;

                }
            }
        }
    }


    fprintf(f, "max error bipectrum continuity = %lf\n", max_error);
    fprintf(f, "max abs of bispecturm = %lf\n\n", max_global);

    max_error = 0;
    max_global = 0;
    for (n1_global = 1; n1_global <= max_check; n1_global = n1_global + 1){
        for (n2_global = 1; n2_global <= max_check; n2_global = n2_global + 1){
            for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
                a = object_checker->check_rotation(10, 5, 2, extended_power_spectrum_for_check, 1.3);
                if (a > max_error)
                    max_error = a;
            }
        }
    }

    fprintf(f, "max error extended power spectrum rotation = %lf\n", max_error);

    for (n1_global = 1; n1_global <= max_check; n1_global = n1_global + 1){
        for (n2_global = 1; n2_global <= max_check; n2_global = n2_global + 1){
            for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
                a = object_checker->check_continuity(10, 3, epsilon_continuity, 1.3, extended_power_spectrum_for_check);
                if (a > max_error)
                    max_error = a;
            }
        }
    }

    fprintf(f, "max error extended power spectrum continuity = %lf\n", max_error);
    fprintf(f, "max abs of extended power spectrum = %lf\n\n", max_global);

    max_error = 0;
    max_global = 0;
    max_re_global = 0;
    max_im_global = 0;

    for (n1_global = 1; n1_global <= max_check; n1_global = n1_global + 1){
        for ( n2_global = 1; n2_global <= max_check; n2_global = n2_global + 1){
           for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
               for (l1_global = 0; l1_global <= max_check; l1_global = l1_global + 1){
                   for (l2_global = 0; l2_global <= max_check; l2_global = l2_global + 1){
                       for (type_global = 0; type_global <= 1; type_global = type_global + 1){
                             a = object_checker->check_rotation(5, 3, 1, extended_bispectrum_for_check, 1.3);
                             if (a > max_error)
                                 max_error = a;
                       }
                   }
               }
           }
        }
    }
    fprintf(f, "max error extended bispectrum rotation = %lf\n", max_error);
    max_error = 0;
    for (n1_global = 1; n1_global <= max_check; n1_global = n1_global + 1){
        for ( n2_global = 1; n2_global <= max_check; n2_global = n2_global + 1){
           for (l_global = 0; l_global <= max_check; l_global = l_global + 1){
               for (l1_global = 0; l1_global <= max_check; l1_global = l1_global + 1){
                   for (l2_global = 0; l2_global <= max_check; l2_global = l2_global + 1){
                      for (type_global = 0; type_global <= 1; type_global = type_global + 1){
                         a = object_checker->check_continuity(10, 1, epsilon_continuity, 1.3, extended_bispectrum_for_check);
                         if (a > max_error)
                              max_error = a;
                       }
                   }

               }
           }
        }
    }
    fprintf(f, "max error extended bispectrum continuity = %lf\n", max_error);
    fprintf(f, "max abs of extended bispecrum real part = %lf\n", max_re_global);
    fprintf(f, "max abs of extended bispecrum imaginary part = %lf\n\n", max_im_global);

    max_error = 0;
    max_global = 0;
    for (n_global = 0; n_global <= max_check; n_global = n_global + 1){
        for (n1_global = 0; n1_global <= max_check; n1_global = n1_global + 1){
            for (n2_global = 0; n2_global <= max_check; n2_global = n2_global + 1){
                a =  object_checker->check_rotation(5, 3, 1, bispectrum4_for_check, 1.3);
                if (a > max_error)
                   max_error = a;
            }
        }
    }

    fprintf(f, "max error bispectrum4 rotation = %lf\n", max_error);

    max_error = 0;
    for (n_global = 0; n_global <= max_check; n_global = n_global + 1){
        for (n1_global = 0; n1_global <= max_check; n1_global = n1_global + 1){
            for (n2_global = 0; n2_global <= max_check; n2_global = n2_global + 1){
                a =  object_checker->check_continuity(10, 1, epsilon_continuity, 1.3, bispectrum4_for_check);
                if (a > max_error)
                   max_error = a;
            }
        }
    }
    fprintf(f, "max error bispectrum4 continuity = %lf\n", max_error);
    fprintf(f, "max abs of bispectrum4 = %lf\n\n", max_global);

    max_error = 0;
    max_global = 0;

    for (n_global = 0; n_global <= max_check; n_global = n_global + 1){
        a = object_checker->check_rotation(5, 3, 1, power_spectrum4_for_check, 1.3);
        if (a > max_error)
            max_error = a;
    }


    fprintf(f, "max error power spectrum 4 rotation = %lf\n", max_error);

    for (n_global = 0; n_global <= max_check; n_global = n_global + 1){
        a = object_checker->check_continuity(10, 1, epsilon_continuity, 1.3, power_spectrum4_for_check);
        if (a > max_error)
            max_error = a;
    }

    fprintf(f, "max error power spectrum 4 continuity = %lf\n", max_error);
    fprintf(f, "max abs of power spectrum 4 = %lf\n\n", max_global);




    max_error = 0;

    for (n_global = 0; n_global <= max_check; n_global = n_global + 1){
        for (n1_global = 0; n1_global <= max_check; n1_global = n1_global + 1){
            for (n2_global = 0; n2_global <= max_check; n2_global = n2_global + 1){
                a =  object_checker->check_compress(5, 3, bispectrum4_for_check, 1.3);
                if (a > max_error)
                   max_error = a;
            }
        }
    }

    fprintf(f, "max error bispectrum4 compress = %lf\n", max_error);

    printf("finished");
    fclose(f);
}



check_implementor::check_implementor()
{
    object_checker = new checker();
}
