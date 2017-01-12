#ifndef HPC_MATVEC_SCAL_H
#define HPC_MATVEC_SCAL_H 1

#include <cassert>
#include <type_traits>
#include <hpc/ulmblas/gescal.h>
#include <hpc/matvec/isgematrix.h>
#include <hpc/matvec/isdensevector.hpp>

namespace hpc { namespace matvec {

template <typename Alpha, typename MA>
typename std::enable_if<IsGeMatrix<MA>::value,
         void>::type
scal(const Alpha &alpha, MA &A)
{
    typedef typename MA::Index  Index;

    const Index m = A.numRows;
    const Index n = A.numCols;

    ulmblas::gescal(m, n, alpha,
                    A.data, A.incRow, A.incCol);
}

template <typename Alpha, typename VX>
typename std::enable_if<IsDenseVector<VX>::value,
         void>::type
scal(const Alpha &alpha, VX &x)
{
    ulmblas::scal(x.length, alpha, x.data, x.inc);
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_SCAL_H
