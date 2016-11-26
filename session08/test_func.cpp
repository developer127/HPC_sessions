#include "test_func.hpp"

namespace test {

double
asumDiffGeMatrix(const GeMatrix& A, const GeMatrix& B)
{
    assert(A.m == B.m && A.n == B.n);

    double asum = 0;

    for (std::size_t j=0; j<A.n; ++j) {
        for (std::size_t i=0; i<A.m; ++i) {
            asum += fabs(B(i,j) - A(i,j));
        }
    }
    return asum;
}

void
printGeMatrixInMemory(GeMatrix& A)
{
    std::size_t length = A.m * A.n;

    for(size_t i= 0; i<length; ++i) {
        fmt::printf("%5.1lf", A.data[i]);
    }
    fmt::printf("\n");
}

}       //test
