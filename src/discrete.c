/*
 * discrete.c
 *
 *  Created on: 2017/06/04
 */
#include "metropolis.h"
#include <math.h>

double multinomialLogPmf(unsigned int *sample_x, double *alpha, int n_dimentions) {
    int i;
    double r = 0.0;
    for(i=0; i < n_dimentions; i++) {
        r += sample_x[i] * log(alpha[i]);
    }
    return r;
}
