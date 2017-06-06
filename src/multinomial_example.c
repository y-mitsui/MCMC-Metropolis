/*
 * multinomial_example.c
 *
 *  Created on: 2017/06/04
 */
#include "metropolis.h"
#include <math.h>
#include <string.h>
#define N_SAMPLES 10000
#define N_DIMENTIONS 100
#define N_ITER 500000

typedef struct {
    unsigned int *sample;
    int n_sample;
    int n_dimentions;
}MultinomialArgs;

double multinomialLogLikelihood(MultinomialArgs *multinomial_args, Parameter *parameters) {
    int i;
    double r = 0.;

    for(i=0; i < multinomial_args->n_dimentions; i++) {
        r += multinomial_args->sample[i] * log(parameters[0].normalized_parameters[i]);
    }
    return r;
}

int main(void) {
    int i, j;
    double true_prob[N_DIMENTIONS];
    double dirichlet_alpha[N_DIMENTIONS];
    unsigned int multi_variate[N_DIMENTIONS];
    gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
    unsigned int *sample = dMalloc(sizeof(unsigned int) * N_DIMENTIONS);
    memset(sample, 0, sizeof(unsigned int) * N_DIMENTIONS);

    for(j=0; j < N_DIMENTIONS; j++) {
        dirichlet_alpha[j] = 1.0;
    }
    gsl_ran_dirichlet(rng, N_DIMENTIONS, dirichlet_alpha, true_prob);

    for(i=0; i < N_SAMPLES; i++) {
        gsl_ran_multinomial(rng, N_DIMENTIONS, 1, true_prob, multi_variate);
        for(j=0; j < N_DIMENTIONS; j++) {
            if(multi_variate[j] == 1)
                sample[j]++;
        }
    }

    printf("true prob [");
	for(j=0; j < N_DIMENTIONS; j++) {
		printf("%f, ", (double)sample[j] / N_SAMPLES);
	}
	printf("]\n");

    MultinomialArgs multinomial_args;

    multinomial_args.sample = sample;
    multinomial_args.n_sample = N_SAMPLES;
    multinomial_args.n_dimentions = N_DIMENTIONS;

    Parameter normal_parameters[1];
    normal_parameters[0].type = SIMPLEX;
    normal_parameters[0].number = N_DIMENTIONS;
    normal_parameters[0].random_scale = 5e-4;
    normal_parameters[0].lower = 0.0;
    normal_parameters[0]._parameters = dMalloc(sizeof(double) * normal_parameters[0].number);
    normal_parameters[0].normalized_parameters = normal_parameters[0]._parameters;

    metropolis((double (*)(void *, Parameter *))multinomialLogLikelihood, &multinomial_args, normal_parameters, sizeof(normal_parameters) / sizeof(normal_parameters[0]), N_ITER, N_ITER / 2, 10000);
    for(int i=0; i < N_DIMENTIONS; i++) {
        printf("%f ", normal_parameters[0].mean_parameters[i]);
    }
    puts("");
    return 0;
}
