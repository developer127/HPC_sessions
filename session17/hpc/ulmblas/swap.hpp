#ifndef HPC_ULMBLAS_SWAP_HPP
#define HPC_ULMBLAS_SWAP_HPP 1

#include <utility>

namespace hpc { namespace ulmblas {

template <typename Index, typename TX, typename TY>
void
swap(Index n, TX *x, Index incX, TY *y, Index incY)
{
    TY tmp;
    for (Index j=0; j<n; ++j) {
        tmp = y[j*incY];
        y[j*incY] = x[j*incX];
        x[j*incX] = tmp;
    }
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_SWAP_HPP
