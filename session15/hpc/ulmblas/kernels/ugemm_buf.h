#ifndef HPC_ULMBLAS_KERNELS_UGEMM_BUF_H
#define HPC_ULMBLAS_KERNELS_UGEMM_BUF_H 1

#include <hpc/ulmblas/blocksize.h>
#include <hpc/ulmblas/geaxpy.h>
#include <hpc/ulmblas/gescal.h>

namespace hpc { namespace ulmblas {

template <typename Index, typename T, typename Beta, typename TC>
void
ugemm(Index kc, T alpha,
      const T *A, const T *B,
      Beta beta,
      TC *C, Index incRowC, Index incColC)
{
    const Index MR = BlockSize<T>::MR;
    const Index NR = BlockSize<T>::NR;
    T P[BlockSize<T>::MR*BlockSize<T>::NR];

    ugemm(kc, alpha, A, B, T(0), P, Index(1), MR);
    gescal(MR, NR, beta, C, incRowC, incColC);
    geaxpy(MR, NR, T(1), P, Index(1), MR, C, incRowC, incColC);
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_KERNELS_UGEMM_BUF_H
