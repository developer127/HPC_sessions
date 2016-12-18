#ifndef HPC_ULMLAPACK_SWAP_HPP
#define HPC_ULMLAPACK_SWAP_HPP 1

#include <cassert>
#include <hpc/ulmblas/swap.hpp>

namespace hpc { namespace ulmlapack {

template <typename Index, typename TA, typename TP>
void
swap(Index m, Index n,
     TA *A, Index incRowA, Index incColA,
     Index k0, Index k1,
     TP *p, Index incP)
{
    assert(0<=k0);
    assert(k0<=k1);
    assert(k1<=m);

    for (Index k = k0; k< k1; ++k) {
        if (p[k*incP]!= k) {
            hpc::ulmblas::swap(n, &A[k*incRowA], incColA,
                                  &A[p[k*incP]*incRowA], incColA);
        }
    }
}

} } // namespace ulmlapack, hpc

#endif // HPC_ULMLAPACK_SWAP_HPP
