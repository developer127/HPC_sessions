#include <cstdlib>
#include <iostream>
#include <hpc/cuda/check.h>
#include <hpc/cuda/properties.h>

#define N 257

__global__ void axpy(std::size_t n, double alpha, double* x, double* y)
{
    std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < n) {
        y[tid] += alpha * x[tid];
    }
}

template <typename T, typename ALPHA>
__global__ void scal(std::size_t n, ALPHA alpha, T* x)
{
    std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < n) {
        x[tid] *= alpha;
    }
}

std::size_t ceildiv(std::size_t x, std::size_t y) {
   /* note that we expect x > 0 && y > 0;
      not safe against overflows but we expect y to be small */
   return (x + y - 1) / y;
}

int main() {
    double a[N];
    for (std::size_t i = 0; i < N; ++i) {
        a[i] = i;
    }

    /* transfer vectors to GPU memory */
    double* cuda_a;
    CHECK_CUDA(cudaMalloc, (void**)&cuda_a, N * sizeof(double));
    CHECK_CUDA(cudaMemcpy, cuda_a, a, N * sizeof(double),
               cudaMemcpyHostToDevice);

    /* execute kernel function on GPU */
    std::size_t warp_size = hpc::cuda::get_warp_size(); /* typically 32 */
    std::size_t nof_warps = ceildiv(N, warp_size);
    std::size_t warps_per_block =
        hpc::cuda::get_max_threads_per_block() / warp_size / 4; /* typically 8 */
    std::size_t nof_blocks = ceildiv(nof_warps, warps_per_block);
    std::size_t threads_per_block;

    if (nof_blocks == 1) {
        threads_per_block = N;
    } else {
        threads_per_block = warps_per_block * warp_size;
    }
    //axpy<<<nof_blocks, threads_per_block>>>(N, 2.0, cuda_a, cuda_b);
    scal<<<nof_blocks, threads_per_block>>>(N,2, cuda_a);

    /* transfer result vector from GPU to host memory */
    CHECK_CUDA(cudaMemcpy, a, cuda_a, N * sizeof(double),
               cudaMemcpyDeviceToHost);
    /* free space allocated at GPU memory */
    CHECK_CUDA(cudaFree, cuda_a);

    /* print result */
    for (std::size_t i = 0; i < N; ++i) {
        std::cout << " " << a[i];
        if (i % 10 == 0) std::cout << std::endl;
    }
    std::cout << std::endl;
}
