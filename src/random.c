/*
 * random.c
 *
 *  Created on: 2017/06/03
 *      Author: yosuke
 */
#include "metropolis.h"
#include <math.h>

double xor128(){
    static unsigned int x=123456789,y=362436069,z=521288629,w=88675123;
    unsigned int t;
    t=(x^(x<<11));
    x=y;
    y=z;
    z=w;
    w=(w^(w>>19))^(t^(t>>8));
    return (double)w/(double)0xFFFFFFFF;
}

double gsl_rng_uniform_pos2(){
    double r;
    do{
        r=xor128();
    }while(r==0.0);
    return r;
}

void randn(double *result, int n){
    double x, y, r2;
    int i;
    for (i=0; i < n; i++) {
        do{
            /* choose x,y in uniform square (-1,-1) to (+1,+1) */
            x = -1 + 2 * gsl_rng_uniform_pos2 ();
            y = -1 + 2 * gsl_rng_uniform_pos2 ();

            /* see if it is in the unit circle */
            r2 = x * x + y * y;
        }while (r2 > 1.0 || r2 == 0);

        /* Box-Muller transform */
        result[i]=y * sqrt (-2.0 * log (r2) / r2);
    }
}
