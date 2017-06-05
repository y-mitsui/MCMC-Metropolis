/*
 * multinomial_example.c
 *
 *  Created on: 2017/06/04
 */
#include "metropolis.h"
#define N_SAMPLES 500
#define N_DIMENTIONS 3
#define N_ITER 1000

typedef struct {
    unsigned int *sample;
    int n_sample;
    int n_dimentions;
}MultinomialArgs;

double multinomialLogLikelihood(MultinomialArgs *multinomial_args, Parameter *parameters) {
    int i;
    double r = 0.;

    for(i=0; i < multinomial_args->n_sample; i++) {
        r += multinomialLogPmf(&multinomial_args->sample[i * multinomial_args->n_dimentions], parameters[0].normalized_parameters, multinomial_args->n_dimentions);
    }
    return r;
}

int main(void) {
    int i, j;
    double true_prob[N_DIMENTIONS] = {0.1, 0.5, 0.4};
    unsigned int multi_variate[N_DIMENTIONS];
    gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
    unsigned int *sample = dMalloc(sizeof(unsigned int) * N_SAMPLES * N_DIMENTIONS);

    for(i=0; i < N_SAMPLES; i++) {
        gsl_ran_multinomial(rng, 3, 1, true_prob, multi_variate);
        for(j=0; j < N_DIMENTIONS; j++) {
            sample[i * N_DIMENTIONS + j] = multi_variate[j];
        }
    }
    MultinomialArgs multinomial_args;

    multinomial_args.sample = sample;
    multinomial_args.n_sample = N_SAMPLES;
    multinomial_args.n_dimentions = N_DIMENTIONS;

    Parameter normal_parameters[1];
    normal_parameters[0].type = SIMPLEX;
    normal_parameters[0].number = N_DIMENTIONS;
    normal_parameters[0].random_scale = 1e-2;
    normal_parameters[0].lower = 0.0;
    normal_parameters[0]._parameters = dMalloc(sizeof(double) * normal_parameters[0].number);
    normal_parameters[0].normalized_parameters = normal_parameters[0]._parameters;


    metropolis((double (*)(void *, Parameter *))multinomialLogLikelihood, &multinomial_args, normal_parameters, sizeof(normal_parameters) / sizeof(normal_parameters[0]), 1000000, 500000, 10000);
    for(int i=0; i < N_DIMENTIONS; i++) {
        printf("%f ", normal_parameters[0].mean_parameters[i]);
    }
    puts("");
    return 0;
}
