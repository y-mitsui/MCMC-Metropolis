/*
 * metropolis.h
 *
 *  Created on: 2017/06/03
 *      Author: yosuke
 */

#ifndef SRC_METROPOLIS_H_
#define SRC_METROPOLIS_H_
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

enum PARAMETER_TYPE {
    SIMPLEX,
    REAL
};

typedef struct {
    enum PARAMETER_TYPE type;
    int number;
    double random_scale;
    double lower;
    double upper;
    double *_parameters;
    double *_adopt_parameters;
    double *normalized_parameters;
    double *mean_parameters;
}Parameter;

void randn(double *result, int n);
double xor128();
void *dMalloc(int size);
double logNormalPdf(gsl_vector *sample_x, double *means, double *covars, int n_dimentions);
double normal1DLogPdf(double sample_x, double mean, double std);
double multinomialLogPmf(unsigned int *sample_x, double *alpha, int n_dimentions);
void metropolis(double (*logLikelihood)(void *, Parameter*), void* args, Parameter* parameters, int n_parameters, int n_iter, int warmup, int print_log_freq);
gsl_matrix *gsl_matrix_clone(const gsl_matrix *src);
gsl_vector *gsl_vector_clone(const gsl_vector *src);
double gsl_det(gsl_matrix *m);


#endif /* SRC_METROPOLIS_H_ */
