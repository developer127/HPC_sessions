#include <stddef.h>
#include <math.h>

double
asumDiff(size_t m, size_t n,
         const double *A,
         size_t incRowA, size_t incColA,
         const double *B,
         size_t incRowB, size_t incColB)
{
   double diff = 0;
   for (size_t i = 0; i < m; ++i) {
      for (size_t j = 0; j < n; ++j) {
         diff += fabs(B[i*incRowB+j*incColB] - A[i*incRowA+j*incColA]);
      }
   }
   return diff;
}

