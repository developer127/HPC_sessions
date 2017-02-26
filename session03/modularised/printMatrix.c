#include <stdio.h>
#include <stddef.h>

void
printMatrix(size_t m, size_t n,
            const double *A,
            size_t incRowA, size_t incColA)
{
   for (size_t i = 0; i < m; ++i) {
      printf("   ");
      for (size_t j = 0; j < n; ++j) {
         printf("%4.1lf ", A[i*incRowA+j*incColA]);
      }
      printf("\n");
   }
   printf("\n");
}
