#ifndef HPC_MATVEC_IAMAX_HPP
#define HPC_MATVEC_IAMAX_HPP 1

#include <cassert>
#include <type_traits>
#include <hpc/ulmblas/iamax.hpp>
#include <hpc/matvec/isdensevector.hpp>

namespace hpc { namespace matvec {

template <typename VX>
typename std::enable_if<IsDenseVector<VX>::value,
         typename VX::Size>::type
iamax(const VX &x)
{
    return ulmblas::iamax(x.length, x.data, x.inc);
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_IAMAX_HPP
