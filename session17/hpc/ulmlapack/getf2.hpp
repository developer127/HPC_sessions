#ifndef HPC_ULMLAPACK_GETF2_H
#define HPC_ULMLAPACK_GETF2_H 1

#include <hpc/ulmblas/iamax.hpp>
#include <hpc/ulmblas/swap.hpp>
#include <hpc/ulmblas/scal.hpp>
#include <hpc/ulmblas/ger.hpp>

namespace hpc { namespace ulmlapack {

template <typename Index, typename TA, typename TP>
std::ptrdiff_t
getf2(Index m, Index n, TA *A, Index incRowA, Index incColA, TP *p, Index incP)
{
    Index i;
    Index mn = (m<n)? m:n;
    for (Index j=0; j<mn; ++j) {
        i = j + hpc::ulmblas::iamax(mn-j, &A[j*incRowA+j*incColA], incRowA);
        if (i!=j) {
            hpc::ulmblas::swap(n, &A[j*incRowA], incColA,
                               &A[i*incRowA], incColA);
        }
        p[j * incP]=i;
        if (A[j*(incRowA+incColA)] == 0) {
            return j;
        }
        hpc::ulmblas::scal(m-j-1,
                           TA(1)/A[j*(incRowA+incColA)],
                           &A[(j+1)*incRowA+j*incColA], incRowA);
        hpc::ulmblas::ger(m-j-1, n-j-1, TA(-1),
                          &A[(j+1)*incRowA+j*incColA], incRowA,
                          &A[j*incRowA+(j+1)*incColA], incColA,
                          &A[(j+1)*(incRowA+incColA)], incRowA, incColA);
    }
    return -1;
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMLAPACK_GETF2_H
