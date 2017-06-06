/*
 * mixture_bernoulli.c
 *
 *  Created on: 2017/06/05
 */
#include "metropolis.h"
#include <stdio.h>
#include <math.h>
#define N_ITER 1000

typedef struct {
    unsigned int *sample;
    int n_sample;
    int n_dimentions;
    int n_components;
}MixtureBernoulliArgs;

double mixtureBernoulliLogLikelihood(MixtureBernoulliArgs *multinomial_args, Parameter *parameters) {
	int i, j, k;
	double r = 0.;

	double *normalized_parameters1 = malloc(sizeof(double) * multinomial_args->n_components * multinomial_args->n_dimentions);
	double *normalized_parameters0 = malloc(sizeof(double) * multinomial_args->n_components * multinomial_args->n_dimentions);
	double *log_weights = malloc(sizeof(double) * multinomial_args->n_components);
	for(j=0; j < multinomial_args->n_components; j++) {
		for(k=0; k < multinomial_args->n_dimentions; k++) {
			normalized_parameters1[j * multinomial_args->n_dimentions + k] = log(parameters[j + multinomial_args->n_sample + 1].normalized_parameters[k]);
			normalized_parameters0[j * multinomial_args->n_dimentions + k] = log(1 - parameters[j + multinomial_args->n_sample + 1].normalized_parameters[k]);
		}
		log_weights[j] = log(parameters[0].normalized_parameters[j]);
	}
	for(i=0; i < multinomial_args->n_sample; i++) {
		for(j=0; j < multinomial_args->n_components; j++) {
			double bernoulli_prob = 0.;
			for(k=0; k < multinomial_args->n_dimentions; k++) {
				int x = multinomial_args->sample[i * multinomial_args->n_dimentions + k];
				bernoulli_prob += (x==1) ? normalized_parameters1[j * multinomial_args->n_dimentions + k] : normalized_parameters0[j * multinomial_args->n_dimentions + k];
			}

			r += parameters[i + 1].normalized_parameters[j] * (log_weights[j] + bernoulli_prob);
		}
	}

	free(normalized_parameters1);
	free(normalized_parameters0);
	free(log_weights);
	return r;
}
double mixtureBernoulliLogLikelihood2(MixtureBernoulliArgs *multinomial_args, Parameter *parameters) {
    int i, j, k;
    double r = 0.;
    double *normalized_parameters1 = malloc(sizeof(double) * multinomial_args->n_components * multinomial_args->n_dimentions);
    double *normalized_parameters0 = malloc(sizeof(double) * multinomial_args->n_components * multinomial_args->n_dimentions);
    double *log_weights = malloc(sizeof(double) * multinomial_args->n_components);
    for(j=0; j < multinomial_args->n_components; j++) {
    	for(k=0; k < multinomial_args->n_dimentions; k++) {
    		normalized_parameters1[j * multinomial_args->n_dimentions + k] = log(parameters[j + 1].normalized_parameters[k]);
    		normalized_parameters0[j * multinomial_args->n_dimentions + k] = log(1 - parameters[j + 1].normalized_parameters[k]);
    	}
    	log_weights[j] = log(parameters[0].normalized_parameters[j]);
    }
    for(i=0; i < multinomial_args->n_sample; i++) {
    	double mixture_bernoulli_prob = 0.;
    	for(j=0; j < multinomial_args->n_components; j++) {
    		double bernoulli_prob = 0;
			for(k=0; k < multinomial_args->n_dimentions; k++) {
				int x = multinomial_args->sample[i * multinomial_args->n_dimentions + k];
				bernoulli_prob += (x==1) ? normalized_parameters1[j * multinomial_args->n_dimentions + k] : normalized_parameters0[j * multinomial_args->n_dimentions + k];
			}
			//mixture_bernoulli_prob += log(parameters[0].normalized_parameters[j]) + bernoulli_prob;
			mixture_bernoulli_prob += log_weights[j] + bernoulli_prob;
    	}
    	//mixture_bernoulli_prob = (mixture_bernoulli_prob > 1e-15) ? mixture_bernoulli_prob : 1e-15;
    	//r += log(mixture_bernoulli_prob);
    	r += mixture_bernoulli_prob;
    }
    free(normalized_parameters1);
    free(normalized_parameters0);
    free(log_weights);
    return r;
}

int mixtureBernoull(unsigned int *sample_X, int n_samples, int n_dimentions, int n_components, double **means) {
	MixtureBernoulliArgs mixture_bernoulli_args;

	mixture_bernoulli_args.sample = sample_X;
	mixture_bernoulli_args.n_sample = n_samples;
	mixture_bernoulli_args.n_dimentions = n_dimentions;
	mixture_bernoulli_args.n_components = n_components;

	printf("n_components + n_samples + 1:%d\n", n_components + n_samples + 1);
	Parameter *normal_parameters = dMalloc(sizeof(Parameter) * (n_components + n_samples + 1));
	normal_parameters[0].type = SIMPLEX;
	normal_parameters[0].number = n_components;
	normal_parameters[0].random_scale = 1e-3;
	normal_parameters[0].lower = 1e-7;
	normal_parameters[0]._parameters = dMalloc(sizeof(double) * normal_parameters[0].number);
	normal_parameters[0].normalized_parameters = normal_parameters[0]._parameters;

	for (int i=1; i <= n_samples + 1; i++) {
		normal_parameters[i].type = ONE_HOT;
		normal_parameters[i].number = n_components;
		normal_parameters[i].random_scale = 5e-3;
		normal_parameters[i].lower = 0.;
		normal_parameters[i].upper = 1.;
		normal_parameters[i]._parameters = dMalloc(sizeof(double) * normal_parameters[i].number);
		normal_parameters[i].normalized_parameters = normal_parameters[i]._parameters;
	}

	for(int i=n_samples + 1; i <= n_samples + n_components + 1; i++) {
		normal_parameters[i].type = REAL;
		normal_parameters[i].number = n_dimentions;
		normal_parameters[i].random_scale = 5e-3;
		normal_parameters[i].lower = 1e-8;
		normal_parameters[i].upper = 1.0 - normal_parameters[i].lower;
		normal_parameters[i]._parameters = dMalloc(sizeof(double) * normal_parameters[i].number);
		printf("1:%p\n", normal_parameters[i]._parameters);
		normal_parameters[i].normalized_parameters = normal_parameters[i]._parameters;
	}
	metropolis((double (*)(void *, Parameter *))mixtureBernoulliLogLikelihood, &mixture_bernoulli_args, normal_parameters, n_components + n_samples + 1, 100000, 50000, 10000);

	for(int i=0; i < n_components; i++) {
		for(int j=0; j < n_dimentions; j++) {
			means[i][j] = normal_parameters[i + n_samples + 1].mean_parameters[j];
		}
		puts("");
	}
	return 0;
}

int main(void) {
	unsigned int sample_X[] = {0,0,1,1,0,0,0,1,0};
	double **mean = malloc(sizeof(double*) * 3);
	int n_dims = 3;
	int n_samples = 3;
	for(int i=0; i  < 3; i++) {
		mean[i] = malloc(sizeof(double) * 3);
	}
	mixtureBernoull(sample_X, n_samples, n_dims, 2, mean);
	return 0;
}

