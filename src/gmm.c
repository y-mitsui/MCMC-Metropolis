/*
 * gmm.c
 *
 *  Created on: 2017/06/04
 */

#include "metropolis.h"
#include <math.h>


double normal1DLogPdf(double sample_x, double mean, double std) {
    return log(1.0 / sqrt(2 * M_PI * std * std)) + -((sample_x - mean) * (sample_x - mean)) / (2 * std * std);
}
double logNormalPdf(gsl_vector *sample_x, double *means, double *covars, int n_dimentions) {
    int i, j;
    gsl_matrix *sigma=gsl_matrix_alloc(n_dimentions, n_dimentions);

    gsl_vector *u=gsl_vector_alloc(n_dimentions);
    for(i=0;i < n_dimentions;i++){
        gsl_vector_set(u, i, means[i]);
    }
    for(i=0;i < n_dimentions;i++){
        gsl_matrix_set(sigma,i, i,fabs(covars[i]));
    }

    double *pp=&covars[n_dimentions];
    for(i = 0;i < n_dimentions;i++){
        for(j = i+1; j < n_dimentions; j++){
            double aa=(fabs(*pp) < 1e-10) ? 1e-10:fabs(*pp);
            gsl_matrix_set(sigma,i,j,aa);
            gsl_matrix_set(sigma,j,i,aa);
            pp++;
        }
    }

    /* 1/det(sigma) */
    double sigmaDetInv=1.0/sqrt(gsl_det(sigma));
    double constant=1.0/pow(sqrt(2*M_PI),(double)n_dimentions);

    /*tmp=inv(sigma)*/
    gsl_matrix *tmp=gsl_matrix_clone(sigma);

    gsl_linalg_cholesky_decomp(tmp);
    gsl_linalg_cholesky_invert(tmp);

    gsl_vector *vecX = gsl_vector_clone(sample_x);
    gsl_vector *vecTmp = gsl_vector_alloc(n_dimentions);

    gsl_vector_sub(vecX,u);
    gsl_blas_dgemv (CblasTrans, 1.0, tmp, vecX,0.0,vecTmp);
    double num;
    gsl_blas_ddot (vecTmp, vecX,&num);
    double r=log(constant*sigmaDetInv*exp(-0.5*num));
    gsl_matrix_free(tmp);
    gsl_matrix_free(sigma);
    gsl_vector_free(u);
    gsl_vector_free(vecX);
    gsl_vector_free(vecTmp);
    return r;
}

double GmmPdf(void *arg,double *parameter){
    /*normalPDFArg *ctx=arg;
    double sum=0.0;
    int i;
    for(i=0;i<ctx->numCorpoment;i++){
        sum+=ctx->n*log(parameter[i])+multiNormLog(ctx,&parameter[i*ctx->numParameter+ctx->numCorpoment]);
    }
    return sum;*/
    return 0;
}
