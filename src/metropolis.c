/*
 * metropolis.c
 *
 *  Created on: 2017/06/03
 */
#include <stdio.h>
#include "metropolis.h"
#define N_ITER 1000

double dmax(double x, double y) {
    return (x > y) ? x : y;
}

void metropolis(double (*logLikelihood)(void *, Parameter*), void* args, Parameter* parameters, int n_parameters, int n_iter, int print_log_freq) {
    int iter, j, k;
    double current_loglikelyfood;

    for(j=0; j < n_parameters; j++) {
        parameters[j]._adopt_parameters = dMalloc(sizeof(double) * parameters[j].number);
        for(k = 0; k < parameters[j].number; k++) {
            parameters[j]._adopt_parameters[k] = 10.;
        }
    }

    for(iter=0; iter < n_iter; iter++) {
        for (j = 0; j < n_parameters; j++) {
            if(parameters[j].type == SIMPLEX) {
                double parameter_sum = 0.0;
                for(k = 0; k < parameters[j].number; k++) {
                    double random_number;
                    randn(&random_number, 1);
                    parameters[j]._parameters[k] = dmax(parameters[j]._adopt_parameters[k] + random_number * parameters[j].random_scale, 1e-5);
                    parameter_sum += parameters[j]._parameters[k];
                }

                for(k = 0; k < parameters[j].number; k++) {
                    parameters[j].normalized_parameters[k] = parameters[j]._parameters[k] / parameter_sum;
                }
            }else if (parameters[j].type == REAL) {
                for(k = 0; k < parameters[j].number; k++) {
                    double random_number;
                    randn(&random_number, 1);
                    parameters[j].normalized_parameters[k] = dmax(parameters[j]._adopt_parameters[k] + random_number * parameters[j].random_scale, parameters[j].lower);
                }
            }
        }

        double condinate_loglikelyfood = logLikelihood(args, parameters);
        const double ratio = exp(condinate_loglikelyfood - current_loglikelyfood);
        if (iter==0 || xor128() < ratio) {
            for (j = 0; j < n_parameters; j++) {
                for(k = 0; k < parameters[j].number; k++) {
                    parameters[j]._adopt_parameters[k] = parameters[j].normalized_parameters[k];
                }
            }
            current_loglikelyfood = condinate_loglikelyfood;
        }

        if(iter % print_log_freq == 0) {
            printf("%d / %d loglikelyfood:%f\n", iter, n_iter, condinate_loglikelyfood);
        }
    }
}

/*
int main2(void) {
    double current_loglikelyfood;
    int iter, j, k;
    int n_parameter_types = 2;
    Parameter gmm_parameters[2];
    gmm_parameters[0].type = SIMPLEX;
    gmm_parameters[0].number = 2;
    gmm_parameters[0]._parameters = dMalloc(sizeof(double) * gmm_parameters[0].number);
    gmm_parameters[0].normalized_parameters = dMalloc(sizeof(double) * gmm_parameters[0].number);
    gmm_parameters[1].type = REAL;
    gmm_parameters[1].number = 2;
    gmm_parameters[1]._parameters = dMalloc(sizeof(double) * gmm_parameters[1].number);
    gmm_parameters[1].normalized_parameters = gmm_parameters[1]._parameters;

    for(iter=0; iter < N_ITER; iter++) {
        for (j = 0; j < n_parameter_types; j++) {
            if(gmm_parameters[j].type == SIMPLEX) {
                double parameter_sum = 0.0;
                for(k = 0; k < gmm_parameters[j].number; k++) {
                    double random_number;
                    randn(&random_number, 1);
                    gmm_parameters[j]._parameters[k] = dmax(random_number + gmm_parameters[j]._adopt_parameters[k], 0);
                    parameter_sum += gmm_parameters[j]._parameters[k];
                }

                for(k = 0; k < gmm_parameters[j].number; k++) {
                    gmm_parameters[j].normalized_parameters[k] = gmm_parameters[j]._parameters[k] / parameter_sum;
                }
            }else if (gmm_parameters[j].type == REAL) {
                for(k = 0; k < gmm_parameters[j].number; k++) {
                    double random_number;
                    randn(&random_number, 1);
                    gmm_parameters[j].normalized_parameters[j] = gmm_parameters[j].normalized_parameters[k] + random_number;
                }
            }
        }

        double condinate_loglikelyfood = gmm_loglikelyfood(gmm_parameters);
        const double ratio = exp(condinate_loglikelyfood - current_loglikelyfood);
        if (iter==0 || xor128() < ratio) {
            for (j = 0; j < n_parameter_types; j++) {
                for(k = 0; k < gmm_parameters[j].number; k++) {
                    gmm_parameters[j]._adopt_parameters[k] = gmm_parameters[j].normalized_parameters[k];
                }
            }
            current_loglikelyfood = condinate_loglikelyfood;
        }
    }
}*/
