#ifndef HPC_ULMLAPACK_GETRF_HPP
#define HPC_ULMLAPACK_GETRF_HPP 1

#include <hpc/ulmlapack/getf2.hpp>

namespace hpc { namespace ulmlapack {

template <typename Index, typename TA, typename TP>
std::ptrdiff_t
getrf(Index m, Index n, TA *A, Index incRowA, Index incColA, TP *p, Index incP)
{
    Index mn   = std::min(m,n);
    Index bs   = 64;
    std::ptrdiff_t info = -1;

    if (bs<=1 || bs>mn) {
	/* matrix too small -> call unblocked variant */
        hpc::ulmlapck::getf2(m, n, A, incRowA, incColA, p, incP);
    } else {
	/* blocked LU factorization */
    }
    return info;
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMLAPACK_GETRF_HPP
