/*
 * common.c
 *
 *  Created on: 2017/06/04
 *      Author: yosuke
 */
#include <stdio.h>
#include <stdlib.h>


void *dMalloc(int size) {
    void *r = malloc(size);
    if (r == NULL) {
        fprintf(stderr, "out of memory\n");
        exit(1);
    }
    return r;
}
