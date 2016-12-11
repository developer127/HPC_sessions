#ifndef HPC_MATVEC_APPLY_H
#define HPC_MATVEC_APPLY_H 1

#include <type_traits>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace matvec {

template <typename MA, typename Func>
typename std::enable_if<IsGeMatrix<MA>::value,
         void>::type
apply(MA &A, Func func)
{
    typedef typename MA::Index    Index;

    if (A.incRow<A.incCol) {
        for (Index j=0; j<A.numCols; ++j) {
            for (Index i=0; i<A.numRows; ++i) {
                func(A(i,j), i, j);
            }
        }
    } else {
        for (Index i=0; i<A.numRows; ++i) {
            for (Index j=0; j<A.numCols; ++j) {
                func(A(i,j), i, j);
            }
        }
    }
}

template <typename MA, typename Func>
typename std::enable_if<IsGeMatrix<MA>::value,
         void>::type
apply(const MA &A, Func func)
{
    typedef typename MA::Index    Index;

    if (A.incRow<A.incCol) {
        for (Index j=0; j<A.numCols; ++j) {
            for (Index i=0; i<A.numRows; ++i) {
                func(A(i,j), i, j);
            }
        }
    } else {
        for (Index i=0; i<A.numRows; ++i) {
            for (Index j=0; j<A.numCols; ++j) {
                func(A(i,j), i, j);
            }
        }
    }
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_APPLY_H
