#ifndef INC_PACK_HPP
#define INC_PACK_HPP 1

#include <hpc/ulmblas/blockSize.hpp>

namespace hpc { namespace ulmblas {

template <typename TA, typename TP, typename Size, typename Index>
void
pack_A(Size mc, Size kc,
       const TA *A, Index incRowA, Index incColA,
             TP *p)
{
    Size MR = BlockSize<TP>::MR;
    Size mp = (mc+MR-1) / MR;

    for (Size j=0; j<kc; ++j) {
        for (Size l=0; l<mp; ++l) {
            for (Size i0=0; i0<MR; ++i0) {
                Size i  = l*MR + i0;
                Size nu = l*MR*kc + j*MR + i0;
                p[nu]   = (i<mc) ? A[i*incRowA+j*incColA]
                                 : TP(0);
            }
        }
    }
}

template <typename TB, typename TP, typename Size, typename Index>
void
pack_B(Size kc, Size nc,
       const TB *B, Index incRowB, Index incColB,
             TP *p)
{
    Size NR = BlockSize<TP>::NR;
    Size np = (nc+NR-1) / NR;

    for (Size l=0; l<np; ++l) {
        for (Size j0=0; j0<NR; ++j0) {
            for (Size i=0; i<kc; ++i) {
                Size j  = l*NR+j0;
                Size nu = l*NR*kc + i*NR + j0;
                p[nu]   = (j<nc) ? B[i*incRowB+j*incColB]
                                 : TP(0);
            }
        }
    }
}

}}          // namespace hpc::ulmblas
#endif      // INC_PACK_HPP
