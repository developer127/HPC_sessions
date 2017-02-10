#ifndef HPC_CUDA_AXPY_HPP
#define HPC_CUDA_AXPY_HPP 1

#include<cassert>

namespace hpc { namespace cuda {

template <typename Index, typename ALPHA, typename TX, typename TY>
__global__ void axpy(Index n, ALPHA alpha, TX* x, TY* y)
{
    assert(n>=0);
    std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < std::size_t(n)) {
        y[tid] += alpha * x[tid];
    }
}

} }     // namespace hpc::cuda
#endif  //HPC_CUDA_AXPY_HPP
