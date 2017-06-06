/*
 * metropolis.c
 *
 *  Created on: 2017/06/03
 */

#include <stdio.h>
#include <string.h>
#include "metropolis.h"
#define N_ITER 1000

double dmax(double x, double y) {
    return (x > y) ? x : y;
}

double dmin(double x, double y) {
    return (x < y) ? x : y;
}

double range(double val, double min, double max) {
    return dmin(dmax(val, min), max);
}

void metropolis(double (*logLikelihood)(void *, Parameter*), void* args, Parameter* parameters, int n_parameters, int n_iter, int warmup, int print_log_freq) {
    int iter, j, k;
    double current_loglikelyfood;

    for(j=0; j < n_parameters; j++) {
        parameters[j]._adopt_parameters = dMalloc(sizeof(double) * parameters[j].number);
        parameters[j].mean_parameters = dMalloc(sizeof(double) * parameters[j].number);
        memset(parameters[j].mean_parameters, 0, sizeof(double) * parameters[j].number);
        for(k = 0; k < parameters[j].number; k++) {
            parameters[j]._adopt_parameters[k] = range(xor128(), parameters[j].lower, parameters[j].upper);
        }
    }
    gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
    for(iter=0; iter < n_iter; iter++) {
        for (j = 0; j < n_parameters; j++) {
            if(parameters[j].type == SIMPLEX) {
                double parameter_sum = 0.0;
                for(k = 0; k < parameters[j].number; k++) {
                    double random_number;
                    randn(&random_number, 1);
                    parameters[j]._parameters[k] = dmax(random_number * parameters[j].random_scale + parameters[j]._adopt_parameters[k] , 1e-8);
                    parameter_sum += parameters[j]._parameters[k];
                }

                for(k = 0; k < parameters[j].number; k++) {
                    parameters[j].normalized_parameters[k] = parameters[j]._parameters[k] / parameter_sum;
                }
            }else if(parameters[j].type == ONE_HOT) {
            	double parameter_sum = 0.0;
				for(k = 0; k < parameters[j].number; k++) {
					double random_number;
					randn(&random_number, 1);
					parameters[j]._parameters[k] = dmax(random_number * parameters[j].random_scale + parameters[j]._adopt_parameters[k] , 1e-8);
					parameter_sum += parameters[j]._parameters[k];
				}

				for(k = 0; k < parameters[j].number; k++) {
					parameters[j]._parameters[k] = parameters[j]._parameters[k] / parameter_sum;
				}
				unsigned int *multi_variate = malloc(sizeof(unsigned int ) * parameters[j].number);
				gsl_ran_multinomial(rng, parameters[j].number, 1, parameters[j]._parameters, multi_variate);
				for(k = 0; k < parameters[j].number; k++) {
					parameters[j].normalized_parameters[k] = (double)multi_variate[k];
				}
				free(multi_variate);

            }else if (parameters[j].type == REAL) {

                for(k = 0; k < parameters[j].number; k++) {
                    double random_number;
                    randn(&random_number, 1);
                    parameters[j].normalized_parameters[k] = range(parameters[j]._adopt_parameters[k] + random_number * parameters[j].random_scale, parameters[j].lower, parameters[j].upper);
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
            printf("%d / %d log likelihood:%f\n", iter, n_iter, condinate_loglikelyfood);
        }
    }
    for (j = 0; j < n_parameters; j++) {
		for(k = 0; k < parameters[j].number; k++) {
			parameters[j].mean_parameters[k] /= iter - warmup;
		}
	}
}
