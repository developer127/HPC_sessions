#ifndef BLAS_H
#define BLAS_H 1

#include <stddef.h>

void gemm_variant1(size_t m, size_t n, size_t k,
                   double alpha,
                   const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                   const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                   double beta,
                   double *C, ptrdiff_t incRowC, ptrdiff_t incColC);

void gemm_variant2(size_t m, size_t n, size_t k,
                   double alpha,
                   const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                   const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                   double beta,
                   double *C, ptrdiff_t incRowC, ptrdiff_t incColC);

void gemm_variant3(size_t m, size_t n, size_t k,
                   double alpha,
                   const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                   const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                   double beta,
                   double *C, ptrdiff_t incRowC, ptrdiff_t incColC);

void gemm_variant4(size_t m, size_t n, size_t k,
                   double alpha,
                   const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                   const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                   double beta,
                   double *C, ptrdiff_t incRowC, ptrdiff_t incColC);

#endif        //BLAS_H
