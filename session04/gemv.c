#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>

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
#define ALPHA 1
#endif

#ifndef BETA
#define BETA 1
#endif

#ifndef T_MIN
#define T_MIN 1
#endif

double A1[MAX_M*MAX_N];
double A2[MAX_M*MAX_N];
double X[MAX_N*INCX];
double Y[MAX_M*INCY];
double Y1[MAX_M*INCY];
double Y2[MAX_M*INCY];

/* return real time in seconds since start of the process */
double
wallTime() {
   static int ticks_per_second = 0;
   if (!ticks_per_second) {
      ticks_per_second = sysconf(_SC_CLK_TCK);
   }
   struct tms timebuf;
   /* times returns the number of real time ticks passed since start */
   return (double) times(&timebuf) / ticks_per_second;
}

void
initMatrix(size_t m, size_t n, double *A, size_t incRowA, size_t incColA) {
   for (size_t j = 0; j < n; ++j) {
      for (size_t i = 0; i < m; ++i) {
         A[i*incRowA+j*incColA] = ((double)rand() - RAND_MAX/2)*200/RAND_MAX;
      }
   }
}

void
copyMatrix(size_t m, size_t n,
         const double *A, size_t incRowA, size_t incColA,
         double *B, size_t incRowB, size_t incColB)
{
   for (size_t j = 0; j < n; ++j) {
      for (size_t i = 0; i < m; ++i) {
         B[i*incRowB+j*incColB] = A[i*incRowA+j*incColA];
      }
   }
}

double
asumDiffMatrix(size_t m, size_t n,
      const double *A, size_t incRowA, size_t incColA,
      const double *B, size_t incRowB, size_t incColB) {
   double diff = 0;
   for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < n; ++j) {
         diff += fabs(B[i*incRowB+j*incColB] - A[i*incRowA+j*incColA]);
      }
   }
   return diff;
}

void
dgemv(size_t m, size_t n,
      double alpha,
      const double *A, size_t incRowA, size_t incColA,
      const double *x, size_t incX,
      double beta,
      double *y, size_t incY) 
{
    size_t i,j;
    if(beta != 1) {
        for(i=0; i<m; ++i){
            y[i*incY] *= beta;
        }

    }
    if(incColA<incRowA){        // case A row major
        for(i=0; i<m; ++i) {
            for(j=0; j<n; ++j) {
                y[i*incY] += alpha * A[i*incRowA+j*incColA] * x[j*incX];
            }
        }
    } else {                    // case A col major
        for(j=0; j<n; ++j) {
            for(i=0; i<m; ++i) {
                y[i*incY] += alpha * A[i*incRowA+j*incColA] * x[j*incX];
            }
        }
    }
}

int
main() {
   initMatrix(MAX_M, MAX_N, A1, 1, MAX_M); /* A1 is in col major */
   initMatrix(MAX_N, 1, X, INCX, 1);
   initMatrix(MAX_M, 1, Y, INCY, 1);

   printf("RUN     M     N");
   printf("  col major  row major");
   printf("  col major  row major");
   printf("      DIFF");
   printf("\n");
   printf("              ");
   printf("   (t in s)   (t in s)");
   printf("   (GFLOPS)   (GFLOPS)");
   printf("         ");
   printf("\n");

   for (size_t i = 0, m = MIN_M, n = MIN_M;
	 m <= MAX_M && n <= MAX_N; ++i, m += 100, n += 100) {
      size_t incRowA, incColA, runs;
      double alpha = ALPHA, beta = BETA;

      /* col major */
      incRowA = 1; incColA = m;
      double t1 = 0;
      runs = 0;
      while (t1 <= T_MIN) {
         copyMatrix(m, 1, Y, INCY, 1, Y1, INCY, 1);
         double t0 = wallTime();
         dgemv(m, n, alpha,
               A1, incRowA, incColA,
               X, INCX,
               beta,
               Y1, INCY);
         t1 += wallTime() - t0;
         ++runs;
      }
      t1 /= runs;

      /* A2 as m x n matrix is identical to A1 as m x n matrix
         but organized in row major: */
      copyMatrix(m, n, A1, 1, m, A2, n, 1);

      /* row major */
      incRowA = n; incColA = 1;
      double t2 = 0;
      runs = 0;
      while (t2 <= T_MIN) {
         copyMatrix(m, 1, Y, INCY, 1, Y2, INCY, 1);
         double t0 = wallTime();
         dgemv(m, n, alpha,
               A2, incRowA, incColA,
               X, INCX,
               beta,
               Y2, INCY);
         t2 += wallTime() - t0;
         ++runs;
      }
      t2 /= runs;

      double diff = asumDiffMatrix(m, 1, Y1, INCY, 1, Y2, INCY, 1);

      printf("%3zu %5zu %5zu ", i, m, n);
      printf("%10.4lf %10.4lf ", t1, t2);
      printf("%10.4lf ", 2*(m/1000.0)*(n/1000.0)/t1);
      printf("%10.4lf ", 2*(m/1000.0)*(n/1000.0)/t2);
      printf("%10.4lf ", diff);
      printf("\n");
   }
}
