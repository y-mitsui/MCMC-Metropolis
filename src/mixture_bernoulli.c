/*
 * mixture_bernoulli.c
 *
 *  Created on: 2017/06/05
 */
#include "metropolis.h"
#define N_SAMPLES 500
#define N_DIMENTIONS 3
#define N_COMPONENTS 2
#define N_ITER 1000

typedef struct {
    unsigned int *sample;
    int n_sample;
    int n_dimentions;
}MixtureBernoulliArgs;

double mixtureBernoulliLogLikelihood(MixtureBernoulliArgs *multinomial_args, Parameter *parameters) {
    int i, j;
    double r = 0.;

    for(i=0; i < multinomial_args->n_sample; i++) {
    	for(j=0; j < multinomial_args->n_dimentions; j++) {
    		int x = multinomial_args->sample[i * multinomial_args->n_dimentions + j];
    		r += x * log(parameters[1].normalized_parameters[j]) + (x - 1) * (1 - parameters[1].normalized_parameters[j]);
    	}
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
	    MixtureBernoulliArgs mixture_bernoulli_args;

	    mixture_bernoulli_args.sample = sample;
	    mixture_bernoulli_args.n_sample = N_SAMPLES;
	    mixture_bernoulli_args.n_dimentions = N_DIMENTIONS;

	    Parameter normal_parameters[1];
	    normal_parameters[0].type = SIMPLEX;
	    normal_parameters[0].number = N_COMPONENTS;
	    normal_parameters[0].random_scale = 1e-2;
	    normal_parameters[0].lower = 0.0;
	    normal_parameters[0]._parameters = dMalloc(sizeof(double) * normal_parameters[0].number);
	    normal_parameters[0].normalized_parameters = normal_parameters[0]._parameters;
	    normal_parameters[1].type = REAL;
		normal_parameters[1].number = N_DIMENTIONS;
		normal_parameters[1].random_scale = 1e-2;
		normal_parameters[1].lower = 0.0;
		normal_parameters[1].upper = 1.0;
		normal_parameters[1]._parameters = dMalloc(sizeof(double) * normal_parameters[0].number);
		normal_parameters[1].normalized_parameters = normal_parameters[0]._parameters;

	    metropolis((double (*)(void *, Parameter *))multinomialLogLikelihood, &mixture_bernoulli_args, normal_parameters, sizeof(normal_parameters) / sizeof(normal_parameters[0]), 100000, 50000, 1000);
	    for(int i=0; i < N_DIMENTIONS; i++) {
	        printf("%f ", normal_parameters[0].mean_parameters[i]);
	    }
	    puts("");
	    return 0;
}
