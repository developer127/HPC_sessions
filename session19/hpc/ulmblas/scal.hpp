#ifndef HPC_ULMBLAS_SCAL_HPP
#define HPC_ULMBLAS_SCAL_HPP 1

namespace hpc { namespace ulmblas {

template <typename Index, typename Alpha, typename TX>
void
scal(Index n,
     const Alpha &alpha,
     TX *x, Index incX)
{
    for (Index j = 0; j<n; ++j) {
        x[j*incX] *= alpha;
    }
}
} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_SCAL_HPP
