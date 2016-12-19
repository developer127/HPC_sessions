#ifndef HPC_ULMBLAS_TRUSM_HPP
#define HPC_ULMBLAS_TRUSM_HPP 1

#include <hpc/ulmblas/geaxpy.h>
#include <hpc/ulmblas/scal.hpp>

namespace hpc { namespace ulmblas {

template <typename Index, typename Alpha, typename TA, typename TB>
void
trusm(Index         m, Index         n,
      const Alpha   &alpha,
      bool          unitDiag,
      const TA      *A, Index incRowA, Index incColA,
      TB            *B, Index incRowB, Index incColB)
{
    if(unitDiag) {
        for (std::ptrdiff_t i = m-2; i>=0; --i) {
            for (std::ptrdiff_t j=m-1; j>i; --j) {
                hpc::ulmblas::geaxpy(Index(1),n,
                       -A[i*incRowA + j*incColA],
                       &B[j*incRowB], incRowB, incColB,
                       &B[i*incRowB], incRowB, incColB);
            }
        }
    } else {
        for (std::ptrdiff_t i = m-1; i>=0; --i) {
            for (std::ptrdiff_t j=m-1; j>i; --j) {
                hpc::ulmblas::geaxpy(Index(1),n,
                       -A[i*incRowA + j*incColA],
                       &B[j*incRowB], incRowB, incColB,
                       &B[i*incRowB], incRowB, incColB);
            }
            hpc::ulmblas::scal(n,
                               TA(1)/A[i*(incRowA+incColA)],
                               &B[i*incRowB], incColB);
        }

    }
}


} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_TRUSM_HPP
