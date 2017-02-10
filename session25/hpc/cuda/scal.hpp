#ifndef HPC_CUDA_SCAL_CUH
#define HPC_CUDA_SCAL_CUH 1

#include<cassert>

namespace hpc { namespace cuda {

template <typename Index, typename ALPHA, typename T>
__global__ void scal(Index n, ALPHA alpha, T* x)
{
    assert(n>=0);
    std::size_t tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < std::size_t(n)) {
        x[tid] *= alpha;
    }
}

} }     // namespace hpc::cuda
#endif  //HPC_CUDA_SCAL_CUH
