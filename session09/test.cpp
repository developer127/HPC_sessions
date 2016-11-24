#include <cstdio>
#include <random>
#include "bench.hpp"
#include <fmt/printf.hpp>
#include "ulmblas.hpp"
#include "gemm_refcolmajor.hpp"
#include <type_traits>

template <typename T, typename Size, typename Index>
void
printGeMatrix(Size m, Size n, const T *A, Index incRowA, Index incColA)
{
    for (Size i=0; i<m; ++i) {
        for (Size j=0; j<n; ++j) {
            fmt::printf(" %f", A[i*incRowA+j*incColA]);
        }
        fmt::printf("\n");
    }
    fmt::printf("\n");
}

#ifndef REAL
#define REAL  double
#endif


template <typename T>
void
checkType(T)
{
    fmt::printf("Unknown type\n");
}

void
checkType(int)
{
    fmt::printf("int\n");
}

void
checkType(double)
{
    fmt::printf("double\n");
}

template <typename T, typename S>
typename std::common_type<T,S>::type
max(T a, S b)
{
    return (a>b) ? a : b;
}

int
main()
{
    checkType(max(1,   1.0));
    checkType(max(1.0, 1  ));
    checkType(max(1,   1  ));
    checkType(max(1.0, 1.0));

    return 0;
}

