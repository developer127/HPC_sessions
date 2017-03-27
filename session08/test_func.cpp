#include "test_func.hpp"

namespace test {

double
asumDiffMatrix(std::size_t m, std::size_t n,
               const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
               double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB)
{

    double asum = 0;

    for (std::size_t j=0; j<n; ++j) {
        for (std::size_t i=0; i<m; ++i) {
            asum += fabs(B[i*incRowB+j*incColB] - A[i*incRowA+j*incColA]);
        }
    }
    return asum;
}


void
printGeMatrixInMemory(std::size_t m, std::size_t n,
                      const double *A)
{
    for(size_t i= 0; i<m*n; ++i) {
        fmt::printf("%6.1lf", A[i]);
    }
    fmt::printf("\n\n");
}

}       //test
