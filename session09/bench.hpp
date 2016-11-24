#ifndef INC_BENCH_H
#define INC_BENCH_H 1

#include <random>
#include <chrono>

namespace bench {

template <typename T>
struct WallTime
{
    std::chrono::high_resolution_clock::time_point  t0;
    std::chrono::high_resolution_clock::duration    elapsed;

    void
    tic()
    {
        t0 = std::chrono::high_resolution_clock::now();
    }

    T
    toc()
    {
        elapsed = std::chrono::high_resolution_clock::now() - t0;
        return std::chrono::duration<T,std::chrono::seconds::period>
               (elapsed).count();
    }
};

template <typename T, typename Size, typename Index>
void
initGeMatrix(Size m, Size n, T *A, Index incRowA, Index incColA)
{
    std::random_device                       random;
    std::mt19937                             mt(random());
    std::uniform_real_distribution<double>   uniform(-100, 100);

    for (Size j=0; j<n; ++j) {
        for (Size i=0; i<m; ++i) {
            A[i*incRowA+j*incColA] = uniform(mt);
        }
    }
}

template <typename T, typename Size, typename Index>
T
asumDiffGeMatrix(Size m, Size n,
                 const T *A, Index incRowA, Index incColA,
                 T *B, Index incRowB, Index incColB)
{
    T asum = T(0);

    for (Size j=0; j<n; ++j) {
        for (Size i=0; i<m; ++i) {
            asum += fabs(B[i*incRowB+j*incColB] - A[i*incRowA+j*incColA]);
        }
    }
    return asum;
}

}   //bench namespace
#endif
