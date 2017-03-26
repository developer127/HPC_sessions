#include "matrix_func.hpp"

void
initMatrix(std::size_t m, std::size_t n,
           double *A, std::size_t incRowA, std::size_t incColA)
{
    for (std::size_t j=0; j<n; ++j) {
        for (std::size_t i=0; i<m; ++i) {
            A[i*incRowA+j*incColA] = ((double)rand() - RAND_MAX/2)*200/RAND_MAX;
        }
    }
}
