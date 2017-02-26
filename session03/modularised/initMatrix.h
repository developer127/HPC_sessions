#ifndef INIT_MATRIX_H
#define INIT_MATRIX_H 1

#include <stddef.h>

void
initMatrix(size_t m, size_t n,
           double *A,
           size_t incRowA, size_t incColA);

#endif        //INIT_MATRIX_H
