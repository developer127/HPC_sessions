#ifndef HPC_MATVEC_COPY_H
#define HPC_MATVEC_COPY_H 1

#include <cassert>
#include <type_traits>
#include <hpc/ulmblas/gecopy.h>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace matvec {

template <typename MA, typename MB>
typename std::enable_if<IsGeMatrix<MA>::value
                     && IsGeMatrix<MB>::value,
         void>::type
copy(const MA &A, MB &B)
{
    assert(A.numRows==B.numRows);
    assert(A.numCols==B.numCols);

    typedef typename std::common_type<typename MA::Index,
                                      typename MB::Index>::type  Index;

    const Index m = B.numRows;
    const Index n = B.numCols;

    ulmblas::gecopy(m, n,
                    A.data, A.incRow, A.incCol,
                    B.data, B.incRow, B.incCol);
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_COPY_H
