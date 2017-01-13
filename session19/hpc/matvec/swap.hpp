#ifndef HPC_MATVEC_SWAP_HPP
#define HPC_MATVEC_SWAP_HPP 1

#include <cassert>
#include <type_traits>
#include <hpc/ulmblas/swap.hpp>
#include <hpc/matvec/isdensevector.hpp>

namespace hpc { namespace matvec {

template <typename VX, typename VY>
typename std::enable_if<IsDenseVector<VX>::value
                        && IsDenseVector<VY>::value,
         void>::type
swap(VX &x, VY &y)
{
    assert(x.length == y.length);
    ulmblas::swap(x.length, x.data, x.inc, y.data, y.inc);
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_SWAP_HPP
