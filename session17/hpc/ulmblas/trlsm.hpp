#ifndef HPC_ULMBLAS_TRLSM_HPP
#define HPC_ULMBLAS_TRLSM_HPP 1

#include <hpc/ulmblas/geaxpy.h>

namespace hpc { namespace ulmblas {

template <typename Index, typename Alpha, typename TA, typename TB>
void
trlsm(Index         m, Index         n,
      const Alpha   &alpha,
      bool          unitDiag,
      const TA      *A, Index incRowA, Index incColA,
      TB            *B, Index incRowB, Index incColB)
{
    if(unitDiag) {
        for (Index i = 1; i<m; ++i) {
            for (Index j=0; j<i; ++j) {
                hpc::ulmblas::geaxpy(Index(1),n,
                       -A[i*incRowA + j*incColA],
                       &B[j*incRowB], incRowB, incColB,
                       &B[i*incRowB], incRowB, incColB);
            }
        }
    }
}


} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_TRLSM_HPP
