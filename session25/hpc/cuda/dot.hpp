#ifndef HPC_CUDA_DOT_HPP
#define HPC_CUDA_DOT_HPP 1

#include<assert>

namespace hpc { namespace cuda {

template<typename Index, typename TX, typename TY, typename T>
__global__ void dot(Index n,
      const TX* x, Index incX, TY* y, Index incY, T* sums)
{
    std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < n) {
        x[tid] *= y[tid];
    }
}

} }     // namespace hpc::cuda
#endif  //HPC_CUDA_DOT_HPP
