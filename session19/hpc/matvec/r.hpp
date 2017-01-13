#ifndef HPC_MATVEC_R_HPP
#define HPC_MATVEC_R_HPP 1

#include <cassert>
#include <type_traits>
#include <hpc/ulmblas/ger.hpp>
#include <hpc/matvec/isdensevector.hpp>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace matvec {

template <typename Alpha, typename VX, typename VY, typename MA>
typename std::enable_if<IsDenseVector<VX>::value
                    && IsDenseVector<VY>::value
                    && IsGeMatrix<MA>::value,void>::type
r(const Alpha alpha, const VX &x, const VY &y, MA A)
{
    assert(x.length==A.numRows && y.length==A.numCols);
    ulmblas::ger(x.length, y.length, alpha,
                 x.data, x.inc,
                 y.data, y.inc,
                 A.data, A.incRow, A.incCol);
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_R_HPP
