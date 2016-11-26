#ifndef HPC_BENCH_H
#define HPC_BENCH_H 1

#include <chrono>
#include <cmath>
#include <complex>
#include <random>

namespace bench {

template <typename T, typename Index>
T
asumDiffGeMatrix(Index m, Index n,
                 const T *A, Index incRowA, Index incColA,
                 T *B, Index incRowB, Index incColB)
{

    T asum = 0;

    for (Index j=0; j<n; ++j) {
        for (Index i=0; i<m; ++i) {
            asum += std::abs(B[i*incRowB+j*incColB] - A[i*incRowA+j*incColA]);
        }
    }
    return asum;
}

template <typename T, typename Index, typename Initvalue>
void
initGeMatrix(Index m, Index n, T *A, Index incRowA, Index incColA,
             Initvalue initValue)
{
    for (Index j=0; j<n; ++j) {
        for (Index i=0; i<m; ++i) {
            A[i*incRowA+j*incColA] = initValue(i, j, m, n);
        }
    }
}

template <typename T>
struct WallTime
{
    void
    tic()
    {
        t0 = std::chrono::high_resolution_clock::now();
    }

    T
    toc()
    {
        using namespace std::chrono;

        elapsed = high_resolution_clock::now() - t0;
        return duration<T,seconds::period>(elapsed).count();
    }

    std::chrono::high_resolution_clock::time_point t0;
    std::chrono::high_resolution_clock::duration   elapsed;
};


} // namespace bench


#endif // HPC_BENCH_H
