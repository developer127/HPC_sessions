#ifndef HPC_CUDA_ASUM_HPP
#define HPC_CUDA_ASUM_HPP 1

#include<cassert>
#include<hpc/cuda/check.h>

#define N 9192

#define THREADS_PER_BLOCK 256
#define NUM_BLOCKS ((N + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK)

namespace hpc { namespace cuda {

template<typename Index, typename T>
__global__ void asum_step(Index n, const T* x, T* sums) {

    __shared__ T sums_per_block[THREADS_PER_BLOCK];
     std::size_t blockIndex = threadIdx.x;

    // Aggregate the sum in the first entry of sums_per_block
    // assuming THREADS_PER_BLOCK is a power of two
    std::size_t length = blockDim.x;
    while (length > 1) {
        length /= 2;
        if (blockIndex < length) {
            sums_per_block[blockIndex] += sums_per_block[length + blockIndex];
        }
        __syncthreads();
    }
    if (blockIndex == 0) {
        sums[blockIdx.x] = sums_per_block[0];
    }
}

template<typename Index, typename T>
T asum(Index n, const T* x) {
    assert(n>=0);
    Index newlength = n;
    T* cuda_x;
    CHECK_CUDA(cudaMalloc, (void**)&cuda_x, n * sizeof(T));
    CHECK_CUDA(cudaMemcpy, cuda_x, x, n * sizeof(T),
               cudaMemcpyHostToDevice);

    T* cuda_sums;
    CHECK_CUDA(cudaMalloc, (void**)&cuda_sums, NUM_BLOCKS * sizeof(T));

    while(std::size_t(newlength) > blockDim.x) {
        newlength = (newlength + blockDim.x-1)/blockDim.x;
        Index newnumBlocks = (newlength+THREADS_PER_BLOCK -1)/THREADS_PER_BLOCK;
        T* cuda_newsum;
        CHECK_CUDA(cudaMalloc, (void**)&cuda_newsum, newnumBlocks * sizeof(T));

        asum_step<<<newnumBlocks, THREADS_PER_BLOCK>>>
                 (newlength, cuda_sums, cuda_newsum);
        cuda_sums = cuda_newsum;
        CHECK_CUDA(cudaFree, cuda_newsum);
    }
    /* transfer result vector from GPU to host memory */
    double sums[THREADS_PER_BLOCK];
    CHECK_CUDA(cudaMemcpy, sums, cuda_sums, THREADS_PER_BLOCK * sizeof(T),
               cudaMemcpyDeviceToHost);

    /* free space allocated at GPU memory */
    CHECK_CUDA(cudaFree, cuda_newsum);
    CHECK_CUDA(cudaFree, cuda_sums);

    return sums[0];
}

} }     // namespace hpc::cuda
#endif  //HPC_CUDA_ASUM_HPP
