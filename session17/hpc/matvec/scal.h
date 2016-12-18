#ifndef HPC_MATVEC_SCAL_H
#define HPC_MATVEC_SCAL_H 1

#include <cassert>
#include <type_traits>
#include <hpc/ulmblas/gescal.h>
#include <hpc/matvec/isgematrix.h>

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

} } // namespace matvec, hpc

#endif // HPC_MATVEC_SCAL_H
