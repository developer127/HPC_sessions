#ifndef HPC_ULMBLAS_GER_HPP
#define HPC_ULMBLAS_GER_HPP 1

#include <utility>

namespace hpc { namespace ulmblas {

template <typename Index, typename Alpha, typename TX, typename TY, typename TA>
void
ger(Index m, Index n, const Alpha &alpha,
    const TX *x, Index incX,
    const TY *y, Index incY,
    TA *A, Index incRowA, Index incColA)
{
    if (incRowA < incColA) {        // ColMajor
        for (Index j = 0; j<n; ++j) {
            for (Index i=0; i<m; ++i) {
                A[i*incRowA + j*incColA] += alpha * x[i*incX] * y[j*incY];
            }
        }
    } else {                        //RowMajor
        for (Index i=0; i<m; ++i) {
            for (Index j = 0; j<n; ++j) {
                A[i*incRowA + j*incColA] += alpha * x[i*incX] * y[j*incY];
            }
        }
    }
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_GER_HPP
