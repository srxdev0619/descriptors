#include "check_implementor.h"

#define max_length 50
#define max_size_input_data 10000

using namespace std;

double distance_segment(vector a, vector b, vector r){

    if (dot_product(b - a, r - a) < 0)
        return distance(a, r);
    if (dot_product(a - b, r - b) < 0)
        return distance(b, r);
    return abs(vector_product(b - a, r - a)) / abs(b - a);
}

double distance_parallelogram(vector p, vector a1, vector a2, vector r){
    vector n = vector_product(a1, a2);
    matrix A = createMatrix(a1, a2, n);
    A = inverse(A);
    vector b = r - p;
    vector k = mul(A, b);
    double k1 = k.x;
    double k2 = k.y;
    double k3 = k.z;
    if ((k1 >= 0) && (k1 <= 1))
        if ((k2 >= 0) && (k2 <= 1)){
            return abs(n) * abs_double(k3);
        }
    double distances[5];
    distances[1] = distance_segment(p, p + a1, r);
    distances[2] = distance_segment(p, p + a2, r);
    distances[3] = distance_segment(p + a1, p + a1 + a2, r);
    distances[4] = distance_segment(p + a2, p + a1 + a2, r);

    double min = distances[1];
    for (int i = 2; i <= 4; i = i + 1)
        if (distances[i] < min)
            min = distances[i];
    return min;

}

double distance_cell(vector p, vector a1, vector a2, vector a3, vector r){
    double distances[7];
    distances[1] = distance_parallelogram(p, a1, a2, r);
    distances[2] = distance_parallelogram(p, a1, a3, r);
    distances[3] = distance_parallelogram(p, a2, a3, r);
    distances[4] = distance_parallelogram(p + a1, a2, a3, r);
    distances[5] = distance_parallelogram(p + a2, a1, a3, r);
    distances[6] = distance_parallelogram(p + a3, a1, a2, r);

    double min = distances[1];
    for (int i = 2; i <= 6; i = i + 1)
        if (distances[i] < min)
            min = distances[i];
    return min;
}

void clear(bool ***considered){
    for (int i = -max_length; i <= max_length; i = i + 1){
        for (int j = -max_length; j <= max_length; j = j + 1){
            for (int k = -max_length; k <= max_length; k = k + 1){
                considered[i][j][k] = false;
            }
        }
    }
}


void add_atoms(double r_cut, int i0, input_data *data_addition, int l, int k, int m, input_data *data_source, vector a1, vector a2, vector a3){
    a1 = mulV(a1, l);
    a2 = mulV(a2, k);
    a3 = mulV(a3, m);
    vector shift = a1 + a2 + a3;

    vector q;
    for (int i = 1; i <= data_source->num_atoms; i = i + 1){
        if (!((l == 0) && (k == 0) && (m == 0) && (i == i0))){
            q = shift + data_source->vectors[i] - data_source->vectors[i0];
            if (abs(q) <= r_cut){
                data_addition->num_atoms = data_addition->num_atoms + 1;
                data_addition->vectors[data_addition->num_atoms] = q;
                data_addition->type[data_addition->num_atoms] = data_source->type[i];
            }
        }
    }


}
void dfs(int l, int k, int m, double r_cut, int i0, input_data *data_source, input_data*data_ans, vector a1, vector a2, vector a3, bool ***considered);

void dfs_next(int l, int k, int m, double r_cut, int i0, input_data *data_source, input_data*data_ans, vector a1, vector a2, vector a3, bool ***considered){
    if (considered[l][k][m]) return;
    vector shift1 = mulV(a1, l);
    vector shift2 = mulV(a2, k);
    vector shift3 = mulV(a3, m);
    vector shift = shift1 + shift2 + shift3;
    double distance = distance_cell(shift, a1, a2 ,a3, data_source->vectors[i0]);
    if (distance <= r_cut)
        dfs(l, k, m, r_cut, i0, data_source, data_ans, a1, a2, a3, considered);
}

void dfs(int l, int k, int m, double r_cut, int i0, input_data *data_source, input_data*data_ans, vector a1, vector a2, vector a3, bool ***considered){
   add_atoms(r_cut, i0, data_ans, l, k, m, data_source, a1, a2, a3);
   considered[l][k][m] = true;
   dfs_next(l + 1, k, m, r_cut, i0, data_source, data_ans, a1, a2, a3, considered);
   dfs_next(l - 1, k, m, r_cut, i0, data_source, data_ans, a1, a2, a3, considered);

   dfs_next(l, k + 1, m, r_cut, i0, data_source, data_ans, a1, a2, a3, considered);
   dfs_next(l, k - 1, m, r_cut, i0, data_source, data_ans, a1, a2, a3, considered);

   dfs_next(l, k, m + 1, r_cut, i0, data_source, data_ans, a1, a2, a3, considered);
   dfs_next(l, k, m - 1, r_cut, i0, data_source, data_ans, a1, a2, a3, considered);

}

void print_descriptors(input_data data, double r_cut, FILE *output){
    descriptor34* object = new descriptor34(7, 6, 2.5, r_cut);
    object->process_data(data);

    for (int n = 1; n <= 7; n = n + 1){
        for (int l = 0; l <= 6; l = l + 1){
            fprintf(output, "%lf ", object->calculate_power_specturm3(n, l));
        }
    }

    for (int n = 1; n <= 3; n = n + 1){
        for (int l = 0; l <= 2; l = l + 1){
            for (int l1 = 0; l1 <= 2; l1 = l1 + 1){
                for (int l2 = 0; l2 <= 2; l2 = l2 + 1){
                    fprintf(output, "%lf ", object->calculate_bispectrum3(n, l, l1, l2));
                }
            }
        }
    }

    for (double j = 0; j <= 2.6; j = j + 0.5){
        fprintf(output, "%lf ", object->calculate_power_spectrum4(j));
    }

    for (double j = 0; j <= 1.6; j = j + 0.5){
        for (double j1 = 0; j1 <= 1.6; j1 = j1 + 0.5){

                fprintf(output, "%lf ", object->calculate_bispectrum4(j, j1, j1));

        }
    }




    fprintf(output, "\n");
}

int main()
{

    /*matrix A;
    vector a1 = {1, 2, 3};
    vector a2 = {7, 8, 9};
    vector a3 = {2, -1, 3};
    A = createMatrix(a1, a2, a3);
    A.print();
    A = inverse(A);
    printf("\n\n");
    A.print();
    return 0;*/

    printf("poscar file name\n");
    char name_input_poscar[1000];
    scanf("%s", name_input_poscar);
    FILE *input_poscar;
    input_poscar = fopen(name_input_poscar, "r");

    printf("energies file\n");
    char name_input_energies[1000];
    scanf("%s", name_input_energies);
    FILE *input_energies;
    input_energies = fopen(name_input_energies, "r");


    bool ***considered;

    considered = new bool**[2 * max_length + 1];
    considered = considered + max_length;

    for (int i = -max_length; i <= max_length; i = i + 1){
        considered[i] = new bool*[2 * max_length + 1];
        considered[i] = considered[i] + max_length;
        for (int j = -max_length; j <= max_length; j = j + 1){
            considered[i][j] = new bool[2 * max_length + 1];
            considered[i][j] = considered[i][j] + max_length;
        }
    }




    char s[1000];
    int num = 0;

    input_data data_result;
    data_result.type = new int[max_size_input_data + 1];
    data_result.vectors = new vector[max_size_input_data + 1];

    input_data data_now;
    data_now.type = new int[max_size_input_data + 1];
    data_now.vectors = new vector[max_size_input_data + 1];

    double energy_now;
    double r_cut_now = 0.6;

    while (!feof(input_poscar)){
        num = num + 1;
        char name_output[1000];
        sprintf(name_output, "results\\info number %d.txt", num);
        FILE* output;
        output = fopen(name_output, "w");

        double scale;
        vector a1, a2, a3;
        int nZn, nO;
        vector q;

        fgets(s, 1000, input_poscar);


        fgets(s, 1000, input_energies);
        sscanf(s, "%lf", &energy_now);

        fscanf(input_poscar, "%lf\n", &scale);
        fscanf(input_poscar, "%lf %lf %lf\n", &a1.x, &a1.y, &a1.z);
        fscanf(input_poscar, "%lf %lf %lf\n", &a2.x, &a2.y, &a2.z);
        fscanf(input_poscar, "%lf %lf %lf\n", &a3.x, &a3.y, &a3.z);
        fgets(s, 1000, input_poscar);
        if (s[0] != 'Z') printf("not as usual!");
        fscanf(input_poscar, "%d %d\n", &nZn, &nO);
        fgets(s, 1000, input_poscar);
        if (s[0] != 'D') printf("not as usual!");


        data_now.num_atoms = 0;
        for (int i = 1; i <= nZn; i = i + 1){
            fscanf(input_poscar, "%lf %lf %lf\n", &q.x, &q.y, &q.z);
            data_now.num_atoms = data_now.num_atoms + 1;
            data_now.vectors[data_now.num_atoms] = q;
            data_now.type[data_now.num_atoms] = 1;
        }

        for (int i = 1; i <= nO; i = i + 1){
            fscanf(input_poscar, "%lf %lf %lf\n", &q.x, &q.y, &q.z);
            data_now.num_atoms = data_now.num_atoms + 1;
            data_now.vectors[data_now.num_atoms] = q;
            data_now.type[data_now.num_atoms] = 2;
        }

        //fprintf(output, "this configuration consist of %d atoms in unit cell, has energy %lf per unit cell, and chosen cut off radius is %lf\n", data_now.num_atoms, energy_now, r_cut_now);
	fprintf(output,"%lf\n",energy_now);
        for (int i = 1; i <= data_now.num_atoms; i = i + 1){
            clear(considered);
            data_result.num_atoms = 0;
            dfs(0, 0, 0, r_cut_now, i, &data_now, &data_result, a1, a2, a3, considered);
            //fprintf(output, "atom number %d with position %lf %lf %lf and type %d has following %d neibours (in relative cartesian coordinates)\n",
            //        i, data_now.vectors[i].x, data_now.vectors[i].y, data_now.vectors[i].z, data_now.type[i], data_result.num_atoms);
            //for (int j = 1; j <= data_result.num_atoms; j = j + 1){
	    //   fprintf(output, "%d %lf %lf %lf\n", data_result.type[j],
	    //           data_result.vectors[j].x, data_result.vectors[j].y, data_result.vectors[j].z);

            //}
	//fprintf(output, "descriptors for neibours of this atom:\n");
            print_descriptors(data_result, r_cut_now, output);
            fprintf(output, "\n");

        }


        fclose(output);
        printf("%d finished\n", num);
        //break;


    }
    fclose(input_poscar);

    printf("%d", num);
    //return a.exec();
}
