#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "initMatrix.h"
#include "walltime.h"
#include "printMatrix.h"
#include "asumDiff.h"


#ifndef MINDIM_M
#define MINDIM_M    1000
#endif

#ifndef MINDIM_N
#define MINDIM_N    1000
#endif

#ifndef MINDIM_K
#define MINDIM_K    1000
#endif

#ifndef MAXDIM_M
#define MAXDIM_M    8000
#endif

#ifndef MAXDIM_N
#define MAXDIM_N    8000
#endif

#ifndef MAXDIM_K
#define MAXDIM_K    8000
#endif

#ifndef INC_M
#define INC_M   100
#endif

#ifndef INC_N
#define INC_N   100
#endif

#ifndef INC_K
#define INC_K   100
#endif

#ifndef MIN_T
#define MIN_T   1
#endif

double buffer1[MAXDIM_M*MAXDIM_N];
double buffer2[MAXDIM_M*MAXDIM_N];

int
main()
{
   printf("   M    N      t1      t2   t2/t1       diff\n");
   printf("          col-maj row-maj\n");
   printf("============================================\n");

   for (size_t m = MINDIM_M, n = MINDIM_N; m < MAXDIM_M && n < MAXDIM_N;
	    m += INC_M, n += INC_N) {
      size_t runs = 0;
      double t1 = 0;
      do {
         double t0 = walltime();
         // col major
         initMatrix(m, n, buffer1, 1, m);
         t1 += walltime() - t0;
         ++runs;
      } while (t1 < MIN_T);
      t1 /= runs;

      runs = 0;
      double t2 = 0;
      do {
         double t0 = walltime();
         // row major
         initMatrix(m, n, buffer2, n, 1);
         t2 += walltime() - t0;
         ++runs;
      } while (t2 < MIN_T);
      t2 /= runs;

      double diff = asumDiff(m, n, buffer1, 1, m, buffer2, n, 1);

      printf("%4zd %4zd %7.2lf %7.2lf %7.2lf %10.2le\n",
	 m, n, t1, t2, t2/t1, diff);
   }
}
