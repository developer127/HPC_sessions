#ifndef HPC_CUDA_ASUM_HPP
#define HPC_CUDA_ASUM_HPP 1

#include<cassert>
#include<hpc/cuda/check.h>

#define N 8192

#define THREADS_PER_BLOCK 256
#define NUM_BLOCKS ((N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK)

namespace hpc { namespace cuda {

template<typename Index, typename TX, typename TY, typename T>
__global__ void asum(Index n, const TX* x, Index incX, T* sums)
{
    assert(n>=0);
    std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;

    // writing simultaniously into registers. Thus no waiting.
    T res;
    if (tid < std::size_t(n)) {
        res = x[tid * incX] * y[tid * incY];
    } else {
        res = 0;
    }

    /* writing to shared variable
    needs to wait until all blocks are finished */
    __shared__ T sums_per_block[THREADS_PER_BLOCK];
    std::size_t blockIndex = threadIdx.x;
    sums_per_block[blockIndex] = res;

    // Aggregate the sum in the first entry of sums_per_block
    // assuming THREADS_PER_BLOCK is a power of two
    std::size_t length = blockDim.x;
    while (length > 1) {
        __syncthreads();
        length /= 2;
        if (blockIndex < length) {
            sums_per_block[blockIndex] += sums_per_block[length + blockIndex];
        }
    }
    if (blockIndex == 0) {
        sums[blockIdx.x] = sums_per_block[0];
    }
}

} }     // namespace hpc::cuda
#endif  //HPC_CUDA_ASUM_HPP
