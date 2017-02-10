#include <cstdio>
#include <cmath>
#include <hpc/cuda/check.h>
#include <hpc/cuda/dot.hpp>

#ifndef N
#define N 9192
#endif

#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 256
#define NUM_BLOCKS ((N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK)
#endif

// Reference implementation
template<typename Index, typename TX, typename TY>
TX dot(Index n, const TX* x, Index incX, TY* y, Index incY) {

   TX res = 0;
   for (Index index = 0; index < n; ++index) {
      res += x[index * incX] * y[index * incY];
   }
   return res;
}


int main()
{
    using namespace hpc::cuda;
    double a[N]; double b[N];
    for (std::size_t i = 0; i < N; ++i) {
        a[i] = i;
        b[i] = i * i;
    }

    /* transfer vectors to GPU memory */
    double* cuda_a;
    CHECK_CUDA(cudaMalloc, (void**)&cuda_a, N * sizeof(double));
    CHECK_CUDA(cudaMemcpy, cuda_a, a, N * sizeof(double),
               cudaMemcpyHostToDevice);
    double* cuda_b;
    CHECK_CUDA(cudaMalloc, (void**)&cuda_b, N * sizeof(double));
    CHECK_CUDA(cudaMemcpy, cuda_b, b, N * sizeof(double),
               cudaMemcpyHostToDevice);

    double* cuda_sums;
    CHECK_CUDA(cudaMalloc, (void**)&cuda_sums, NUM_BLOCKS * sizeof(double));

    /* execute kernel function on GPU */
    dot<<<NUM_BLOCKS, THREADS_PER_BLOCK>>>(N, cuda_a, 1, cuda_b, 1, cuda_sums);

    /* transfer result vector from GPU to host memory */
    double sums[NUM_BLOCKS];
    CHECK_CUDA(cudaMemcpy, sums, cuda_sums, NUM_BLOCKS * sizeof(double),
               cudaMemcpyDeviceToHost);
    /* free space allocated at GPU memory */
    CHECK_CUDA(cudaFree, cuda_a);
    CHECK_CUDA(cudaFree, cuda_b);
    CHECK_CUDA(cudaFree, cuda_sums);

    /* aggregate block sums */
    double sum = 0;
    for (std::size_t i = 0; i < NUM_BLOCKS; ++i) {
        sum += sums[i];
    }

    /* print difference to local result */
    double local_sum = dot(N, a, 1, b, 1);
    std::printf("diff: %12.4lg\n", std::abs(sum - local_sum));
}
