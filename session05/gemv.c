#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
#include <mkl_types.h>

#ifndef MIN_M
#define MIN_M 1000
#endif

#ifndef MIN_N
#define MIN_N 1000
#endif

#ifndef MAX_M
#define MAX_M 10000
#endif

#ifndef MAX_N
#define MAX_N 10000
#endif

#ifndef INCX
#define INCX 1
#endif

#ifndef INCY
#define INCY 1
#endif

#ifndef ALPHA
#define ALPHA 1.5
#endif

#ifndef BETA
#define BETA 1.5
#endif

#ifndef T_MIN
#define T_MIN 1
#endif

#ifndef DGEMV_DOTF_FUSE
#define DGEMV_DOTF_FUSE  2
#endif

#ifndef DGEMV_AXPYF_FUSE
#define DGEMV_AXPYF_FUSE  2
#endif

double A[MAX_M*MAX_N];
double X[MAX_N*INCX];
double Y[MAX_M*INCY];
double Y1[MAX_M*INCY];
double Y2[MAX_M*INCY];
double Y3[MAX_M*INCY];
double Y4[MAX_M*INCY];
double Y5[MAX_M*INCY];

double
walltime()
{
   struct tms    ts;
   static double ClockTick=0.0;

   if (ClockTick==0.0) {
        ClockTick = 1.0 / ((double) sysconf(_SC_CLK_TCK));
   }
   return ((double) times(&ts)) * ClockTick;
}

void
initMatrix(long m, long n, double *A, long incRowA, long incColA)
{
    int i, j;
    for (j=0; j<n; ++j) {
        for (i=0; i<m; ++i) {
            A[i*incRowA+j*incColA] = ((double)rand() - RAND_MAX/2)*200/RAND_MAX;
        }
    }
}

void
copyMatrix(long m, long n,
           const double *A, long incRowA, long incColA,
           double *B, long incRowB, long incColB)
{
    int i, j;
    for (j=0; j<n; ++j) {
        for (i=0; i<m; ++i) {
            B[i*incRowB+j*incColB] = A[i*incRowA+j*incColA];
        }
    }
}

double
asumDiffMatrix(long m, long n,
               const double *A, long incRowA, long incColA,
               double *B, long incRowB, long incColB)
{
    int i, j;

    double asum = 0;

    for (j=0; j<n; ++j) {
        for (i=0; i<m; ++i) {
            asum += fabs(B[i*incRowB+j*incColB] - A[i*incRowA+j*incColA]);
        }
    }
    return asum;
}

//--------------------BLAS LEVEL1-----------------------------------------------

/** dscal_ulm:  scaling a vector
  * 1.  lengt
  * 2.  beta    scaling factor
  * 3.  *vec    pointer to the first value
  * 4.  incV    increment for the next value of the vector
  *
  * If beta is not zero a well initialized vector is asumed
  * beta = 0 initializes the vector to zero
  */

void
dscal_ulm(size_t length, double beta, double *vec, size_t incV)
{
    if (beta != 1) {
        if(beta != 0) {
            for(size_t i=0; i<length; ++i) {
                vec[i*incV] *= beta;
            }
        } else {
            for(size_t i=0; i<length; ++i) {
                vec[i*incV] = 0;
            }
        }
    }
}


/**
  * ddot_ulm:   vector product u * v
  * 1.  length
  * 2.  *u      ptr to the first vector u
  * 3.  incU    increment to get to the next entry
  * 4.  *v      ptr to the second vector v
  * 5.  incV    increment to get to the next entry
  *
  * we asume that the length of u an v are equal.
  */

double
ddot_ulm(size_t length, const double *u, size_t incU,
         const double *v, size_t incV)
{
    double solution = 0;
    for(size_t i=0; i<length; ++i) {
        solution += u[i*incU] * v[i*incV];
    }
    return solution;
}


/**
  * dotf_ulm:       Matrix_row_panel vektor product
  *                 y <- y + alpha A_panel * x
  * 1.  length      lenght of vector x
  * 2.  alpha       scaling factor
  * 3.  *A_panel    ptr to the first element of the row_panel
  * 4.  incRowA
  * 5.  incColA
  * 6.  *x          ptr to the vector x
  * 7.  incX        increment to get to the next entry
  * 8.  *y          ptr to the solution vector y
  * 9.  incY
  */

void
dotf_ulm(size_t length, double alpha,
         const double *A_panel, size_t incRowA, size_t incColA,
         const double *x, size_t incX,
         double *y, size_t incY)
{
    if(alpha != 0){
        for(size_t j=0; j<length; ++j) {
            for(size_t i=0; i< DGEMV_DOTF_FUSE; ++i) {
                y[i*incY] += alpha * A_panel[i*incRowA+j*incColA] * x[j*incX];
            }
        }
    }
}


/**
  * daxpy_ulm:      y <- y + alpha * x;
  * 1. length
  * 2. alpha        scaling factor
  * 3. *x           ptr to vec. x
  * 4. incX         increment to get to the next element
  * 5. *y           ptr to vec. y
  * 6. incY         increment to get to the next element
  */

void
daxpy_ulm(size_t length, double alpha,
          const double *x, size_t incX,
          double *y, size_t incY)
{
    if(alpha != 0) {
        for(size_t i=0; i<length; ++i) {
            y[i*incY] += alpha * x[i*incX];
        }
    }
}


/**
  * axpyf_ulm:      y <- y + alpha * A_Panel*x;
  * 1. length
  * 2. alpha        scaling factor
  * 3. *A_panel     ptr to the first element of the col_panel
  * 4. incRowA
  * 5. incColA
  * 6. *x           ptr to vec. x
  * 7. incX         increment to get to the next element
  * 8. *y           ptr to vec. y
  * 9. incY         increment to get to the next element
  */

void
axpyf_ulm(size_t length, double alpha,
          const double *A_panel, size_t incRowA, size_t incColA,
          const double *x, size_t incX,
          double *y, size_t incY)
{
    if(alpha != 0) {
        for(size_t i=0; i<length; ++i) {
            double dot = 0;
            for(size_t j=0; j<DGEMV_AXPYF_FUSE; ++j) {
                dot += A_panel[i*incRowA + j*incColA]*x[j*incX];
            }
            y[i*incY] += alpha * dot;
        }
    }
}
//-----------------Blas level2--------------------------------------------------

void
dgemv(const char *trans,
      const MKL_INT *m, const MKL_INT *n, const double *alpha,
      const double *A, const MKL_INT *ldA, const double *x,
      const MKL_INT *incX,
      const double *beta, double *y, const MKL_INT *incY);

void
dgemv_mkl(MKL_INT m, MKL_INT n,
          double alpha,
          const double *A, MKL_INT incRowA, MKL_INT incColA,
          const double *x, MKL_INT incX,
          double beta,
          double *y, MKL_INT incY)
{
    MKL_INT ldA   = (incRowA==1) ? incColA : incRowA;
    char    trans = (incRowA==1) ? 'N' : 'T';
    MKL_INT M     = (incRowA==1) ? m : n;
    MKL_INT N     = (incRowA==1) ? n : m;

    dgemv(&trans, &M, &N, &alpha, A, &ldA, x, &incX, &beta, y, &incY);
}

//------------------------------------------------------------------------------

void
dgemv_ulm(long m, long n,
          double alpha,
          const double *A, long incRowA, long incColA,
          const double *x, long incX,
          double beta,
          double *y, long incY)
{
    long i, j;

    if (beta!=1.0) {
        if(beta != 0.0) {
            for (i=0; i<m; ++i) {
                y[i*incY] *= beta;
            }
        } else {
            for (i=0; i<m; ++i) {
                y[i*incY] = 0;
            }
        }
    }
    if (incRowA > incColA) {
        for (i=0; i<m; ++i) {
            for (j=0; j<n; ++j) {
                y[i*incY] += alpha*A[i*incRowA+j*incColA]*x[j*incX];
            }
        }
    } else  {
        for (j=0; j<n; ++j) {
            for (i=0; i<m; ++i) {
                y[i*incY] += alpha*A[i*incRowA+j*incColA]*x[j*incX];
            }
        }
    }
}

void
dgemv_ulm_mat(size_t m, size_t n,
              double alpha,
              const double *A, size_t incRowA, size_t incColA,
              const double *x, size_t incX,
              double beta,
              double *y, size_t incY)
{
    dscal_ulm(m, beta, y, incY);
    if(incRowA > incColA) {     // A row major
        for(size_t i=0; i<m; ++i) {
            y[i*incY] += alpha * ddot_ulm(n, &A[i*incRowA], incColA, x, incX);
        }
    } else {                    // A col major
        for(size_t j=0; j<n; ++j) {
            daxpy_ulm(m, alpha * x[j*incX],
                      &A[j*incColA], incRowA,
                      y, incY);
        }
    }
}

void
gemvf_ulm(size_t m, size_t n,
          double alpha,
          const double *A, size_t incRowA, size_t incColA,
          const double *x, size_t incX,
          double beta,
          double *y, size_t incY)
{
    dscal_ulm(m, beta, y, incY);

    if(incRowA > incColA) {     // A row major
        size_t mb = m/DGEMV_DOTF_FUSE;
        for(size_t i=0; i<mb; ++i) {
            dotf_ulm(n, alpha,
                     &A[i*incRowA*DGEMV_DOTF_FUSE], incRowA, incColA,
                     x, incX,
                     &y[i*incY*DGEMV_DOTF_FUSE], incY);
        }
        if(n%DGEMV_DOTF_FUSE != 0) {
            dgemv_ulm_mat(m%DGEMV_DOTF_FUSE,n,
                          alpha,
                          &A[(mb)*DGEMV_DOTF_FUSE*incRowA], incRowA, incColA,
                          x, incX, 1,
                          &y[(mb)*DGEMV_DOTF_FUSE*incY], incY);
        }

    } else {                    // A col major
        size_t nb = n/DGEMV_AXPYF_FUSE;
        for(size_t j=0; j<nb; ++j) {
            axpyf_ulm(m, alpha,
                      &A[j*incColA*DGEMV_AXPYF_FUSE], incRowA, incColA,
                      &x[j*incX*DGEMV_AXPYF_FUSE], incX,
                      y, incY);
        }
        if(n%DGEMV_AXPYF_FUSE != 0) {
            dgemv_ulm_mat(m,n%DGEMV_AXPYF_FUSE,
                          alpha,
                          &A[(nb)*DGEMV_AXPYF_FUSE*incColA], incRowA, incColA,
                          &x[(nb)*DGEMV_AXPYF_FUSE*incX], incX, 1, y, incY);
        }
    }
}


//------------------------------------------------------------------------------

#ifndef COLMAJOR
//#define COLMAJOR 1
#define COLMAJOR 0
#endif

int
main()
{
    long runs, i, m, n, incRowA, incColA;
    double t0, t1, t2, t3;
    double diff2;
    double alpha = ALPHA;
    double beta  = BETA;

    initMatrix(MAX_M, MAX_N, A, 1, MAX_M);
    initMatrix(MAX_N, 1, X, INCX, 1);
    initMatrix(MAX_M, 1, Y, INCY, 1);

    printf("# COLMAJOR    = %d\n", COLMAJOR);
    printf("# T_MIN       = %d\n", T_MIN);
    printf("#RUN    M     N  INCROW  INCCOL");
    printf("    GEMV_MKL    GEMV_ULM    GEMVF   ");
    printf("    GEMV_MKL    GEMV_ULM    GEMVF   ");
    printf("       DIFF2");
    printf("\n");
    printf("#                              ");
    printf("    (t in s)    (t in s)    (t in s)");
    printf("    (MFLOPS)    (MFLOPS)    (MFLOPS)");
    printf("           ");
    printf("\n");

    for (i=0, m=MIN_M, n=MIN_N; m<=MAX_M && n<=MAX_N; ++i, m+=100, n+=100) {

        if (COLMAJOR) {
            incRowA = 1;
            incColA = m;
        } else {
            incRowA = n;
            incColA = 1;
        }

        t1   = 0;
        runs = 0;
        do {
            copyMatrix(MAX_M, 1, Y, INCY, 1, Y1, INCY, 1);
            t0 = walltime();
            dgemv_mkl(m, n, alpha,
                      A, incRowA, incColA,
                      X, INCX,
                      beta,
                      Y1, INCY);
            t1 += walltime() - t0;
            ++runs;
        } while (t1<T_MIN);
        t1 /= runs;

        t2   = 0;
        runs = 0;
        do {
            copyMatrix(MAX_M, 1, Y, INCY, 1, Y2, INCY, 1);
            t0 = walltime();
            dgemv_ulm(m, n, alpha,
                      A, incRowA, incColA,
                      X, INCX,
                      beta,
                      Y2, INCY);
            t2 += walltime() - t0;
            ++runs;
        } while (t2<T_MIN);
        t2 /= runs;

        t3   = 0;
        runs = 0;
        do {
            copyMatrix(MAX_M, 1, Y, INCY, 1, Y3, INCY, 1);
            t0 = walltime();
            gemvf_ulm(m, n, alpha,
                      A, incRowA, incColA,
                      X, INCX,
                      beta,
                      Y3, INCY);
            t3 += walltime() - t0;
            ++runs;
        } while (t3<T_MIN);
        t3 /= runs;


        diff2 = asumDiffMatrix(m, 1, Y1, INCY, 1, Y3, INCY, 1);

        printf("%3ld %5ld %5ld %7ld %7ld ", i, m, n, incRowA, incColA);
        printf("%11.4lf %11.4lf %11.4lf" , t1, t2, t3);
        printf("%11.4lf ", 2*(m/1000.0)*(n/1000.0)/t1);
        printf("%11.4lf ", 2*(m/1000.0)*(n/1000.0)/t2);
        printf("%11.4lf ", 2*(m/1000.0)*(n/1000.0)/t3);
        printf("%11.4lf ", diff2);
        printf("\n");
    }

    return 0;
}
