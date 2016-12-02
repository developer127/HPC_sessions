#ifndef HPC_ULMBLAS_MGEMM_H
#define HPC_ULMBLAS_MGEMM_H 1

#include <algorithm>
#include <hpc/ulmblas/blocksize.h>
#include <hpc/ulmblas/geaxpy.h>
#include <hpc/ulmblas/gescal.h>
#include <hpc/ulmblas/kernels/ugemm.h>

namespace hpc { namespace ulmblas {

template <typename Index, typename T, typename Beta, typename TC>
void
mgemm(Index mc, Index nc, Index kc,
      T alpha,
      const T *A, const T *B,
      Beta beta,
      TC *C, Index incRowC, Index incColC)
{
    const Index MR = BlockSize<T>::MR;
    const Index NR = BlockSize<T>::NR;
    T C_[BlockSize<T>::MR*BlockSize<T>::NR];

    const Index mp  = (mc+MR-1) / MR;
    const Index np  = (nc+NR-1) / NR;
    const Index mr_ = mc % MR;
    const Index nr_ = nc % NR;

    for (Index j=0; j<np; ++j) {
        const Index nr = (j!=np-1 || nr_==0) ? NR : nr_;

        for (Index i=0; i<mp; ++i) {
            const Index mr = (i!=mp-1 || mr_==0) ? MR : mr_;

            if (mr==MR && nr==NR) {
                ugemm(kc, alpha,
                      &A[i*kc*MR], &B[j*kc*NR],
                      beta,
                      &C[i*MR*incRowC+j*NR*incColC],
                      incRowC, incColC);
            } else {
                std::fill_n(C_, MR*NR, T(0));
                ugemm(kc, alpha,
                      &A[i*kc*MR], &B[j*kc*NR],
                      T(0),
                      C_, Index(1), MR);
                gescal(mr, nr, beta,
                       &C[i*MR*incRowC+j*NR*incColC],
                       incRowC, incColC);
                geaxpy(mr, nr, T(1), C_, Index(1), MR,
                       &C[i*MR*incRowC+j*NR*incColC],
                       incRowC, incColC);
            }
        }
    }
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_MGEMM_H
