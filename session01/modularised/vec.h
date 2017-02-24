#ifndef MODULARISED_VEC_H
#define MODULARISED_VEC_H 1

#include <stddef.h>

double vec_sum(double * vec, size_t lenght, ptrdiff_t inc);
void init_vector(double *vec, size_t len, ptrdiff_t incr);
void print_vector(double *vec, size_t len, ptrdiff_t inc);

#endif      //MODULARISED_VEC_H
