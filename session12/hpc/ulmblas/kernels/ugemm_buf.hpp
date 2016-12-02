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

    for (Size l=0; l<MR*NR; ++l) {
        P[l] = T(0);
    }
    for (Size l=0; l<kc; ++l) {
        for (Size j=0; j<NR; ++j) {
            for (Size i=0; i<MR; ++i) {
                P[i+j*MR] += A[i+l*MR]*B[l*NR+j];
            }
        }
    }
    gescal(MR, NR, beta, C, incRowC, incColC);
    for (Size j=0; j<NR; ++j) {
        for (Size i=0; i<MR; ++i) {
            C[i*incRowC+j*incColC] += alpha*P[i+j*MR];
        }
    }
}

}}}     // namespace hpc::ulmblas::kernels
#endif  // INC_ULMBLAS_KERNELS_UGEMM_BUF_HPP
