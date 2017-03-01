#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
#include <math.h>
#include "blas.h"
#include "initMatrix.h"
#include "walltime.h"
#include "printMatrix.h"
#include "asumDiff.h"

#ifndef MINDIM_M
#define MINDIM_M    50
#endif

#ifndef MINDIM_N
#define MINDIM_N    50
#endif

#ifndef MINDIM_K
#define MINDIM_K    50
#endif

#ifndef MAXDIM_M
#define MAXDIM_M    2000
#endif

#ifndef MAXDIM_N
#define MAXDIM_N    2000
#endif

#ifndef MAXDIM_K
#define MAXDIM_K    2000
#endif

#ifndef INC_M
#define INC_M   50
#endif

#ifndef INC_N
#define INC_N   50
#endif

#ifndef INC_K
#define INC_K   50
#endif

#ifndef MIN_T
#define MIN_T   1
#endif

double buffer1[MAXDIM_M*MAXDIM_K];
double buffer2[MAXDIM_K*MAXDIM_N];
double buffer3[MAXDIM_M*MAXDIM_N];
double buffer4[MAXDIM_M*MAXDIM_N];

int
main()
{
    printf("   M    N      t1      t2   t2/t1       diff\n");
    printf("          col-maj row-maj\n");
    printf("============================================\n");

    for (size_t m = MINDIM_M, n = MINDIM_N, k = MINDIM_K;
         m < MAXDIM_M && n < MAXDIM_N && k < MAXDIM_K;
         m += INC_M, n += INC_N, k += INC_K) {
       size_t runs = 0;
       double t1 = 0;
       // col major
       initMatrix(m, n, buffer3, 1, m);     //C
       initMatrix(m, k, buffer1, 1, m);     //A
       initMatrix(k, n, buffer2, 1, k);     //B
       do {
          double t0 = walltime();
          // all row major
          gemm_variant3(m, n, k, 1,
                        buffer1, 1, m,
                        buffer2, 1, k,
                        0,
                        buffer3, 1, m);
          t1 += walltime() - t0;
          ++runs;
       } while (t1 < MIN_T);
       t1 /= runs;

       runs = 0;
       double t2 = 0;
       // col major
       initMatrix(m, n, buffer4, 1, m);     //C
       //initMatrix(m, k, buffer1, 1, m);     //A
       //initMatrix(k, n, buffer2, 1, k);     //B
       do {
          double t0 = walltime();
          // all col major
          gemm_variant4(m, n, k, 1,
                        buffer1, 1, m,
                        buffer2, 1, k,
                        0,
                        buffer4, 1, m);
          t2 += walltime() - t0;
          ++runs;
       } while (t2 < MIN_T);
       t2 /= runs;

       double diff = asumDiff(m, n,
                              buffer3, 1, m,
                              buffer4, 1, m);

       printf("%4zd %4zd %7.2lf %7.2lf %7.2lf %10.2le\n",
              m, n, t1, t2, t2/t1, diff);
    }
}
