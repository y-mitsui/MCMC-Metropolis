/*
 * metropolis.c
 *
 *  Created on: 2017/06/03
 *      Author: yosuke
 */

#include <stdio.h>
#include "metropolis.h"
#define N_ITER 1000

double dmax(double x, double y) {
    return (x > y) ? x : y;
}

void metropolis(double (*logLikelihood)(void *, Parameter*), void* args, Parameter* parameters, int n_parameters, int n_iter, int warmup, int print_log_freq) {
    int iter, j, k;
    double current_loglikelyfood;

    for(j=0; j < n_parameters; j++) {
        parameters[j]._adopt_parameters = dMalloc(sizeof(double) * parameters[j].number);
        parameters[j].mean_parameters = dMalloc(sizeof(double) * parameters[j].number);
        memset(parameters[j].mean_parameters, 0, sizeof(double) * parameters[j].number);
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
                    parameters[j]._parameters[k] = dmax(random_number * parameters[j].random_scale + parameters[j]._adopt_parameters[k] , 1e-5);
                    parameter_sum += parameters[j]._parameters[k];
                }

                for(k = 0; k < parameters[j].number; k++) {
                    parameters[j].normalized_parameters[k] = parameters[j]._parameters[k] / parameter_sum;
                }
            }else if (parameters[j].type == REAL) {
                for(k = 0; k < parameters[j].number; k++) {
                    double random_number;
                    randn(&random_number, 1);
                    parameters[j].normalized_parameters[k] = dmin(dmax(parameters[j]._adopt_parameters[k] + random_number * parameters[j].random_scale, parameters[j].lower), parameters[j].upper);
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

        if(iter > warmup) {
			for (j = 0; j < n_parameters; j++) {
				for(k = 0; k < parameters[j].number; k++) {
					parameters[j].mean_parameters[k] += parameters[j]._adopt_parameters[k];
				}
			}
        }

        if(iter % print_log_freq == 0) {
            printf("%d / %d loglikelyfood:%f\n", iter, n_iter, condinate_loglikelyfood);
        }
    }
    for (j = 0; j < n_parameters; j++) {
		for(k = 0; k < parameters[j].number; k++) {
			parameters[j].mean_parameters[k] /= iter - warmup;
		}
	}
}
