#include <stddef.h>

void
initMatrix(size_t m, size_t n,
           double *A,
           size_t incRowA, size_t incColA)
{
   for (size_t j = 0; j < n; ++j) {
      for (size_t i = 0; i < m; ++i) {
         A[i*incRowA+j*incColA] = (double)j*(double)m+(double)i+1;
      }
   }
}
