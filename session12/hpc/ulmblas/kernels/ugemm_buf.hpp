#ifndef INC_ULMBLAS_KERNELS_UGEMM_BUF_HPP
#define INC_ULMBLAS_KERNELS_UGEMM_BUF_HPP 1

#include <hpc/ulmblas/blockSize.hpp>
#include <hpc/ulmblas/gescal.hpp>

namespace hpc { namespace ulmblas { namespace kernels {

template <typename T, typename TC, typename BETA, typename Size, typename Index>
void
ugemm(Size kc, T alpha,
          const T *A, const T *B,
          BETA beta,
          TC *C, Index incRowC, Index incColC)
{
    const Size MR = BlockSize<T>::MR;
    const Size NR = BlockSize<T>::NR;
    T P[BlockSize<T>::MR*BlockSize<T>::NR];

    ugemm(kc, alpha, A, B, T(0), P, Index(1), Index(MR));
    gescal(MR, NR, beta, C, incRowC, incColC);
    geaxpy(MR, NR, T(1), P, Index(1), Index(MR), C, incRowC, incColC);
}

}}}     // namespace hpc::ulmblas::kernels
#endif  // INC_ULMBLAS_KERNELS_UGEMM_BUF_HPP
