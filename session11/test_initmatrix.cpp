#include <cstdlib>
#include <fmt/printf.hpp>
#include "gematrix.h"
#include <complex>
#include <cmath>

template<typename T, typename Index = std::size_t>
struct AscendingValue {
    Index m, n;

    AscendingValue(Index m, Index n):
        m(m), n(n) {}

    T operator()(Index i, Index j) {
        return T(j)*T(n) + T(i) + T(1);
    }
};

template<typename T, typename Index = std::size_t>
std::complex<T> init_complex(Index i, Index j) {
    return std::complex<T>(T(i),T(j)) + std::complex<T>(T(1),T(1));
}

template<typename T, typename Index = std::size_t>
struct RandomValues {
    std::mt19937 mt;
    std::uniform_real_distribution<T> uniform;

    RandomValues() :
        mt(std::random_device()()), uniform(-100, 100){}

    RandomValues(T seed) :
        mt(seed), uniform(-100, 100){}

    T operator()(Index i, Index j) {
        return uniform(mt);
    }
};

template<typename T, typename Index = std::size_t>
struct IncreasingRandomValues {
    std::mt19937 mt;
    std::uniform_real_distribution<T> uniform;
    T randomNumber;

    IncreasingRandomValues(Index m, Index n) :
        mt(std::random_device()()), uniform(0, 300/(m*n)), randomNumber(-100){}

    T operator()(Index i, Index j) {
        randomNumber += uniform(mt);
        return randomNumber;
    }
};

template<typename Matrix>
void print_matrix(const Matrix& A) {
    for (std::size_t i = 0; i < A.m; ++i) {
        fmt::printf("  ");
        for (std::size_t j = 0; j < A.n; ++j) {
            fmt::printf(" %6.1lf", A(i, j));
        }
        fmt::printf("\n");
    }
}

int main() {
    using namespace matvec;
    GeMatrix<double> A(3, 7, StorageOrder::ColMajor);
    initGeMatrix(A, [&] (std::size_t i, std::size_t j) -> double{
                 return A.m *j + i +1;
    });
    fmt::printf("A:\n");
    print_matrix(A);
    GeMatrix<std::complex<double>> B(3, 5, StorageOrder::ColMajor);
    initGeMatrix(B, init_complex<double>);
    fmt::printf("B:\n");
    print_matrix(B);

    std::random_device random;
    std::mt19937 mt(random());
    std::uniform_real_distribution<double> uniform(-100, 100);
    GeMatrix<double> C(5, 5, StorageOrder::ColMajor);
    initGeMatrix(C, [=] (std::size_t i, std::size_t j) mutable -> double {
        return uniform(mt);
    });

    fmt::printf("C:\n");
    print_matrix(C);
    double maxValue = 0;

    auto maxElement = [&] (std::size_t i, std::size_t j, double ref) -> void {
        if(abs(ref)>maxValue) {
            maxValue = abs(ref);
        }
    };
    fmt::printf("%5.1lf\n",maxValue);
    maxElement(1,2,23);
    fmt::printf("%5.1lf\n",maxValue);
    applyGeMatrix(C, [&] (std::size_t i, std::size_t j, double ref) -> void {
        if(fabs(ref)>maxValue) {
            maxValue = fabs(ref);
        }
    });

    fmt::printf("Max absValue of C is: %6.1lf\n", maxValue);

    GeMatrix<double> D(4, 8, StorageOrder::ColMajor);
    std::uniform_real_distribution<double> uniform_2(0, 300/(C.m*C.n));
    double randomNumber = -100.0;
    initGeMatrix(D, [=,&randomNumber] (std::size_t i, std::size_t j) mutable
                 -> double {
        return randomNumber += uniform_2(mt);
    });
    fmt::printf("D:\n");
    print_matrix(D);

    GeMatrix<double> E(3, 7, StorageOrder::ColMajor);
    applyGeMatrix(E, [&] (std::size_t i, std::size_t j, double ref) -> void {
        E(i,j) = E.m *j + i +1;
    });
    fmt::printf("E:\n");
    print_matrix(E);

    double sumE = 0;
    applyGeMatrix(E, [&] (std::size_t i, std::size_t j, double ref) mutable
                  -> void {
        sumE += ref;
    });
    fmt::printf("\nThe sum of all Elements of E: %6.1lf\n", sumE);

}
