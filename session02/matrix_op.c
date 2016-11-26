#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

#ifndef DIMMAX_M
#define DIMMAX_M 10
#endif
#ifndef DIMMAX_N
#define DIMMAX_N 15
#endif


/**
  * Initialiseation of a matrix in full storage format:
  * 1. m        -> number of rows
  * 2. n        -> number of columns
  * 3. A        -> pointer to the first element
  * 4. incRowA  -> increment to get to the next row
  * 5. incColA  -> increment to get to the next column
  */

void initMatrix(size_t m, size_t n,
                double* A,
                ptrdiff_t incRowA, ptrdiff_t incColA)
{
    for(size_t i=0; i<m; ++i){
        for(size_t j = 0; j<n; ++j){
            A[i*incRowA + j*incColA] = i*n+1+j;
        }
    }
}

void
printMatrix(size_t m, size_t n,
            const double * A, ptrdiff_t incRowA, ptrdiff_t incColA)
{
    for(size_t i=0; i<m; ++i){
        for(size_t j = 0; j<n; ++j){
            printf("%8.2lf ", A[i*incRowA + j*incColA]);
        }
        printf("\n");
    }
}

int main()
{
    size_t m = 7;
    size_t n = 8;

    double* A = (double*) malloc(DIMMAX_M * DIMMAX_N *sizeof(double));
    if(!A){
        printf("allocating memory not successful\n");
        return 1;
    }
    // col major
    printf("Matrix A stored col major: \n");
    initMatrix(m, n,
               A,
               1, m);
    printMatrix(m, n,
                A,
                1, m);
    printf("The transposed matrix A:\n");
    printMatrix(n,m,
                A,
                m,1);
    printf("Content of the memory for A:\n");
    printMatrix(1,n*m,
                A,
                n*m,1);
    printf("Second row of A:\n");
    printMatrix(1,n,
                &A[1],
                1,m);
    printf("Third column of A:\n");
    printMatrix(m,1,
                &A[2*m],
                1,m);
    printf("A block of A:\n");
    printMatrix(2,3,
                &A[3*m+1],
                1,m);


    // row major
    printf("\nMatrix A stored row major: \n");
    initMatrix(m, n,
               A,
               n, 1);
    printMatrix(m, n,
                A,
                n, 1);
    printf("The transposed matrix A:\n");
    printMatrix(n,m,
                A,
                1,n);
    printf("Content of the memory for A:\n");
    printMatrix(1,n*m,
                A,
                n*m,1);
    printf("Second row of A:\n");
    printMatrix(1,n,
                &A[n],
                n,1);
    printf("Third column of A:\n");
    printMatrix(m,1,
                &A[2],
                n,1);
    printf("A block of A:\n");
    printMatrix(2,3,
                &A[1*n+3],
                n,1);


    free(A);
    return 0;
}
