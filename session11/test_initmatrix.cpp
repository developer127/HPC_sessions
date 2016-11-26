#include <cstdlib>
#include <printf.hpp>
#include "gematrix.h"
#include <complex>

template<typename T, typename Index>
T init_value(Index i, Index j, Index m, Index n) {
    return T(j)*T(n) + T(i) + T(1);
}

template<typename T, typename Index>
std::complex<T> init_complex(Index i, Index j, Index m, Index n) {
    return std::complex<T>(T(i),T(j)) + std::complex<T>(T(1),T(1));
}

template<typename Matrix>
void print_matrix(const Matrix& A) {
    for (std::size_t i = 0; i < A.m; ++i) {
        fmt::printf("  ");
        for (std::size_t j = 0; j < A.n; ++j) {
            fmt::printf(" %4.1lf", A(i, j));
        }
        fmt::printf("\n");
    }
}

int main() {
    using namespace matvec;
    GeMatrix<double> A(3, 7, StorageOrder::ColMajor);
    initGeMatrix(A, init_value<double,std::size_t>);
    fmt::printf("A:\n");
    print_matrix(A);
    GeMatrix<std::complex<double>> B(3, 5, StorageOrder::ColMajor);
    initGeMatrix(B, init_complex<double,std::size_t>);
    fmt::printf("B:\n");
    print_matrix(B);
}
