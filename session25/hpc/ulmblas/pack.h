#ifndef HPC_ULMBLAS_PACK_H
#define HPC_ULMBLAS_PACK_H 1

#include <complex>
#include <type_traits>
#include <hpc/ulmblas/blocksize.h>
#include <hpc/ulmblas/pack.h>

namespace hpc { namespace ulmblas {

template <typename Index, typename TA, typename T>
void
pack_A(Index mc, Index kc,
       const TA *A, Index incRowA, Index incColA,
       T *p)
{
    Index MR = BlockSize<T>::MR;
    Index mp = (mc+MR-1) / MR;

    for (Index j=0; j<kc; ++j) {
        for (Index l=0; l<mp; ++l) {
            for (Index i0=0; i0<MR; ++i0) {
                Index i  = l*MR + i0;
                Index nu = l*MR*kc + j*MR + i0;
                p[nu]   = (i<mc) ? A[i*incRowA+j*incColA]
                                 : T(0);
            }
        }
    }
}

template <typename Index, typename TB, typename T>
void
pack_B(Index kc, Index nc,
       const TB *B, Index incRowB, Index incColB,
       T *p)
{
    Index NR = BlockSize<T>::NR;
    Index np = (nc+NR-1) / NR;

    for (Index l=0; l<np; ++l) {
        for (Index j0=0; j0<NR; ++j0) {
            for (Index i=0; i<kc; ++i) {
                Index j  = l*NR+j0;
                Index nu = l*NR*kc + i*NR + j0;
                p[nu]   = (j<nc) ? B[i*incRowB+j*incColB]
                                 : T(0);
            }
        }
    }
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_PACK_H
