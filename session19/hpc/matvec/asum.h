#ifndef HPC_MATVEC_ASUM_H
#define HPC_MATVEC_ASUM_H 1

#include <cassert>
#include <cmath>
#include <type_traits>
#include <hpc/aux/primitive_type.h>
#include <hpc/ulmblas/geaxpy.h>
#include <hpc/matvec/isgematrix.h>

namespace hpc { namespace matvec {

template <typename MA>
typename hpc::aux::PrimitiveType<typename MA::ElementType>::Type
asum(const MA &A)
{
    typedef typename MA::Index                        Index;
    typedef typename MA::ElementType                  T;
    typedef typename hpc::aux::PrimitiveType<T>::Type PT;

    PT result = 0;
    apply(A, [&result](T val, Index, Index)
             {
                 result += std::abs(val);
             });
    return result;
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_ASUM_H
