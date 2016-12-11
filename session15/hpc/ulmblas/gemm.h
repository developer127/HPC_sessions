#ifndef HPC_ULMBLAS_GEMM_H
#define HPC_ULMBLAS_GEMM_H 1

#include <complex>
#include <type_traits>
#include <hpc/aux/isfundamental.h>
#include <hpc/ulmblas/blocksize.h>
#include <hpc/ulmblas/gescal.h>
#include <hpc/ulmblas/pack.h>
#include <hpc/ulmblas/mgemm.h>
#include <hpc/ulmblas/kernels/ugemm.h>
#include <cstdlib>

#ifdef MEM_ALIGN
#   include <xmmintrin.h>
#endif

namespace hpc { namespace ulmblas {

//-----------------------------------------------------------------------------

template <typename T>
typename std::enable_if<hpc::aux::IsFundamental<T>::value,
         T *>::type
malloc(size_t n)
{
#ifdef MEM_ALIGN
    return reinterpret_cast<T *>(_mm_malloc(n*sizeof(T), MEM_ALIGN));
#   else
    return new T[n];
#   endif
}

template <typename T>
typename std::enable_if<! hpc::aux::IsFundamental<T>::value,
         T *>::type
malloc(size_t n)
{
    return new T[n];
}

template <typename T>
typename std::enable_if<hpc::aux::IsFundamental<T>::value,
         void>::type
free(T *block)
{
#ifdef MEM_ALIGN
    _mm_free(reinterpret_cast<void *>(block));
#   else
    delete [] block;
#   endif
}

template <typename T>
typename std::enable_if<! hpc::aux::IsFundamental<T>::value,
         void>::type
free(T *block)
{
    delete [] block;
}

//-----------------------------------------------------------------------------

template <typename Index, typename Alpha,
         typename TA, typename TB,
         typename Beta,
         typename TC>
void
gemm(Index m, Index n, Index k,
     Alpha alpha,
     const TA *A, Index incRowA, Index incColA,
     const TB *B, Index incRowB, Index incColB,
     Beta beta,
     TC *C, Index incRowC, Index incColC)
{
    typedef typename std::common_type<Alpha, TA, TB>::type  T;

    const Index MC = BlockSize<T>::MC;
    const Index NC = BlockSize<T>::NC;
    const Index KC = BlockSize<T>::KC;

    const Index MR = BlockSize<T>::MR;
    const Index NR = BlockSize<T>::NR;

    const Index mb = (m+MC-1) / MC;
    const Index nb = (n+NC-1) / NC;
    const Index kb = (k+KC-1) / KC;

    const Index mc_ = m % MC;
    const Index nc_ = n % NC;
    const Index kc_ = k % KC;

    T *A_ = malloc<T>(MC*KC + MR);
    T *B_ = malloc<T>(KC*NC + NR);

    if (alpha==Alpha(0) || k==0) {
        gescal(m, n, beta, C, incRowC, incColC);
        return;
    }

    for (Index j=0; j<nb; ++j) {
        Index nc = (j!=nb-1 || nc_==0) ? NC : nc_;

        for (Index l=0; l<kb; ++l) {
            Index   kc  = (l!=kb-1 || kc_==0) ? KC : kc_;
            Beta beta_  = (l==0) ? beta : Beta(1);

            pack_B(kc, nc,
                   &B[l*KC*incRowB+j*NC*incColB],
                   incRowB, incColB,
                   B_);

            for (Index i=0; i<mb; ++i) {
                Index mc = (i!=mb-1 || mc_==0) ? MC : mc_;

                pack_A(mc, kc,
                       &A[i*MC*incRowA+l*KC*incColA],
                       incRowA, incColA,
                       A_);

                mgemm(mc, nc, kc,
                      T(alpha), A_, B_, beta_,
                      &C[i*MC*incRowC+j*NC*incColC],
                      incRowC, incColC);
            }
        }
    }
    free(A_);
    free(B_);
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_GEMM_H
