#ifndef INC_MGEMM_HPP
#define INC_MGEMM_HPP 1

#include <hpc/ulmblas/blockSize.hpp>
#include <hpc/ulmblas/kernels/ugemm.hpp>
#include <hpc/ulmblas/gescal.hpp>
#include <hpc/ulmblas/geaxpy.hpp>
#include <hpc/ulmblas/kernels/ugemm.hpp>


namespace hpc { namespace ulmblas {

template <typename T, typename TC, typename BETA, typename Size, typename Index>
void
mgemm(Size mc, Size nc, Size kc,
      T alpha,
      const T *A, const T *B,
      BETA beta,
      TC *C, Index incRowC, Index incColC)
{
    Size MR = BlockSize<TC>::MR;
    Size NR = BlockSize<TC>::NR;
    TC C_[BlockSize<TC>::MR*BlockSize<TC>::NR];

    Size mp  = (mc+MR-1) / MR;
    Size np  = (nc+NR-1) / NR;
    Size mr_ = mc % MR;
    Size nr_ = nc % NR;

    for (Size j=0; j<np; ++j) {
        Size nr = (j!=np-1 || nr_==0) ? NR : nr_;

        for (Size i=0; i<mp; ++i) {
            Size mr = (i!=mp-1 || mr_==0) ? MR : mr_;

            if (mr==MR && nr==NR) {
                kernels::ugemm(kc, alpha,
                               &A[i*kc*MR], &B[j*kc*NR],
                               beta,
                               &C[i*MR*incRowC+j*NR*incColC],
                               incRowC, incColC);
            } else {
                kernels::ugemm(kc, alpha,
                               &A[i*kc*MR], &B[j*kc*NR],
                               TC(0),
                               C_, Index(1), MR);
                gescal(mr, nr, beta,
                       &C[i*MR*incRowC+j*NR*incColC],
                       incRowC, incColC);
                geaxpy(mr, nr, TC(1), C_, Index(1), MR,
                       &C[i*MR*incRowC+j*NR*incColC],
                       incRowC, incColC);
            }
        }
    }
}

}}      // namespace hpc::ulmblas
#endif  // INC_MGEMM_HPP
