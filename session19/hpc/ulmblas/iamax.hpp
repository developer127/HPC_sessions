#ifndef HPC_ULMBLAS_IAMAX_HPP
#define HPC_ULMBLAS_IAMAX_HPP 1

#include <cmath>
#include <fmt/printf.hpp>

namespace hpc { namespace ulmblas {

template <typename Index, typename TX>
Index
iamax(Index n, const TX *x, Index incX)
{
    Index i = 0;
    TX tmp = std::abs(x[0]);            // not save, assume n>0
    for (Index j=1; j<n; ++j) {
        if(std::abs(x[j*incX])>tmp) {
            tmp = std::abs(x[j*incX]);
            i = j;
        }
    }
    return i;
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_IAMAX_HPP
