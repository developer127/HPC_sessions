#ifndef HPC_ULMLAPACK_TF2_HPP
#define HPC_ULMLAPACK_TF2_HPP 1

#include <algorithm>
#include <hpc/matvec/r.hpp>
#include <hpc/matvec/iamax.hpp>
#include <hpc/matvec/scal.h>
#include <hpc/matvec/swap.hpp>
#include <hpc/matvec/isgematrix.h>
#include <hpc/matvec/isdensevector.hpp>
#include <hpc/matvec/densevector.hpp>

namespace hpc { namespace ulmlapack {

template <typename MA, typename VP>
typename std::enable_if<hpc::matvec::IsGeMatrix<MA>::value
                     && hpc::matvec::IsDenseVector<VP>::value,
         typename std::remove_reference<MA>::type::Index>::type
tf2(MA &&A, VP &&p)
{
    typedef typename std::remove_reference<MA>::type MatrixA;

    //typedef typename MatrixA::ElementType    T;
    typedef typename MatrixA::Size           Size;
    //typedef typename MatrixA::Index          Index;

    using namespace hpc::matvec;

    Size m = A.numRows;
    Size n = A.numCols;

    assert(m == p.length);

    Size mn = (m<n)? m:n;

    for(Size j=0; j<mn; ++j) {
        auto aj = A.col(j)(j,m-j);
        p(j) = j + iamax(aj);
        auto ajt = A.row(j);
        if(j!=p(j)) {
            auto apjt = A.row(p(j));
            swap(ajt, apjt);
        }
        if(ajt(j)==0) {
            return j;
        }
        ajt = ajt(j+1,n-j-1);
        scal(1.0/ajt(j), ajt);
        r(1,aj(1,m-j-1),ajt(j+1,n-j-1), A(j+1,j+1,m-j-1,n-j-1));
    }
    return -1;
}

} }     // namespace hpc::ulmlapack
#endif  // HPC_ULMLAPACK_TF2_HPP
