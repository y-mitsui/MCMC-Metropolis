/*
 * normal_example.c
 *
 *  Created on: 2017/06/04
 */
#include "metropolis.h"
#define N_SAMPLES 200
#define N_DIMENTIONS 1
#define N_ITER 1000

typedef struct {
    gsl_vector **sample;
    int n_sample;
}NormalArgs;

double normalLogLikelihood(NormalArgs *normal_args, Parameter *parameters) {
    int i;
    double r = 0.;

    for(i=0; i < normal_args->n_sample; i++) {
        r += normal1DLogPdf(gsl_vector_get(normal_args->sample[i], 0), parameters[0].normalized_parameters[0], parameters[1].normalized_parameters[0]);
    }
    return r;
}

int main(void) {
    int i, j;
    NormalArgs normal_args;
    gsl_vector **sample = malloc(sizeof(gsl_vector*) * N_SAMPLES);

    gsl_rng *rng = gsl_rng_alloc (gsl_rng_default);
    for (i=0; i < N_SAMPLES; i++) {
        sample[i] = gsl_vector_alloc(N_DIMENTIONS);
        for (j=0; j < N_DIMENTIONS; j++) {
            gsl_vector_set(sample[i], j, gsl_ran_gaussian(rng, 10.0) + 10);
        }
    }

    normal_args.sample = sample;
    normal_args.n_sample = N_SAMPLES;

    Parameter normal_parameters[2];
    normal_parameters[0].type = REAL;
    normal_parameters[0].number = 1;
    normal_parameters[0].random_scale = 1e-2;
    normal_parameters[0].lower = -1e+5;
    normal_parameters[0]._parameters = dMalloc(sizeof(double) * normal_parameters[0].number);
    normal_parameters[0].normalized_parameters = normal_parameters[0]._parameters;
    normal_parameters[1].type = REAL;
    normal_parameters[1].number = 1;
    normal_parameters[1].random_scale = 1e-2;
    normal_parameters[1].lower = 1e-2;
    normal_parameters[1]._parameters = dMalloc(sizeof(double) * normal_parameters[1].number);
    normal_parameters[1].normalized_parameters = normal_parameters[1]._parameters;

    metropolis((double (*)(void *, Parameter *))normalLogLikelihood, &normal_args, normal_parameters, sizeof(normal_parameters) / sizeof(normal_parameters[0]), 50000, 25000, 1000);
    printf("%f %f\n", normal_parameters[0].normalized_parameters[0], normal_parameters[1].normalized_parameters[0]);

}
