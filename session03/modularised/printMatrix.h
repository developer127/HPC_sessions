#ifndef PRINT_MATRIX_H
#define PRINT_MATRIX_H 1

#include <stddef.h>

void
printMatrix(size_t m, size_t n,
            const double *A,
            size_t incRowA, size_t incColA);

#endif        //PRINT_MATRIX_H
