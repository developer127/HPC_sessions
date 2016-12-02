#ifndef INC_ULMBLAS_GEMM_HPP
#define INC_ULMBLAS_GEMM_HPP 1

#include <type_traits>
#include <hpc/ulmblas/blockSize.hpp>
#include <hpc/ulmblas/gescal.hpp>
#include <hpc/ulmblas/pack.hpp>
#include <hpc/ulmblas/mgemm.hpp>


namespace hpc { namespace ulmblas {

template <typename TA, typename TB, typename TC, typename ALPHA, typename BETA,
          typename Size, typename Index>
void
gemm(Size m, Size n, Size k,
     ALPHA alpha,
     const TA *A, Index incRowA, Index incColA,
     const TB *B, Index incRowB, Index incColB,
     BETA beta,
     TC *C, Index incRowC, Index incColC)
{
    typedef typename std::common_type<TA,TB,ALPHA>::type TCOMM;

    Size MC = BlockSize<TCOMM>::MC;
    Size NC = BlockSize<TCOMM>::NC;
    Size KC = BlockSize<TCOMM>::KC;

    Size mb = (m+MC-1) / MC;
    Size nb = (n+NC-1) / NC;
    Size kb = (k+KC-1) / KC;

    Size mc_ = m % MC;
    Size nc_ = n % NC;
    Size kc_ = k % KC;

    TCOMM *A_ = new TCOMM[MC*KC];
    TCOMM *B_ = new TCOMM[KC*NC];

    if (alpha==ALPHA(0) || k==0) {
        gescal(m, n, beta, C, incRowC, incColC);
        return;
    }

    for (Size j=0; j<nb; ++j) {
        Size nc = (j!=nb-1 || nc_==0) ? NC : nc_;

        for (Size l=0; l<kb; ++l) {
            Size   kc    = (l!=kb-1 || kc_==0) ? KC : kc_;
            BETA beta_ = (l==0) ? beta : BETA(1);

            pack_B(kc, nc,
                   &B[l*KC*incRowB+j*NC*incColB],
                   incRowB, incColB,
                   B_);

            for (Size i=0; i<mb; ++i) {
                Size mc = (i!=mb-1 || mc_==0) ? MC : mc_;

                pack_A(mc, kc,
                       &A[i*MC*incRowA+l*KC*incColA],
                       incRowA, incColA,
                       A_);

                mgemm(mc, nc, kc,
                      alpha, A_, B_, beta_,
                      &C[i*MC*incRowC+j*NC*incColC],
                      incRowC, incColC);
            }
        }
    }
    delete [] A_;
    delete [] B_;
}

}}      // namespace hpc::ulmblas
#endif  // INC_ULMBLAS_GEMM_HPP
