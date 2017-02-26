#include <cstdio>
#include <cmath>
#include <hpc/cuda/check.h>
#include <hpc/cuda/asum.hpp>

#ifndef N
#define N 9192
#endif

#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 256
#define NUM_BLOCKS ((N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK)
#endif

// Reference implementation
template<typename Index, typename T>
T asum_loc(Index n, const T* x) {
   T res = 0;
   for (Index index = 0; index < n; ++index) {
      res += x[index];
   }
   return res;
}


int main()
{
    using namespace hpc::cuda;
    double a[N];
    for (std::size_t i = 0; i < N; ++i) {
        a[i] = i;
    }

    double sum;

    /* execute kernel function on GPU */
    sum = asum(N, a);

    /* print difference to local result */
    double local_sum = asum_loc(N, a);
    std::printf("diff: %12.4lg\n", std::abs(sum - local_sum));
}
