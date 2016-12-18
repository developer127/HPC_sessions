#ifndef HPC_ULMBLAS_GESCAL_H
#define HPC_ULMBLAS_GESCAL_H 1

namespace hpc { namespace ulmblas {

template <typename Index, typename Alpha, typename TX>
void
gescal(Index m, Index n,
       const Alpha &alpha,
       TX *X, Index incRowX, Index incColX)
{
    if (alpha!=Alpha(1)) {
        if (incRowX<incColX) {
            for (Index j=0; j<n; ++j) {
                for (Index i=0; i<m; ++i) {
                    X[i*incRowX+j*incColX] *= alpha;
                }
            }
        } else {
            for (Index i=0; i<m; ++i) {
                for (Index j=0; j<n; ++j) {
                    X[i*incRowX+j*incColX] *= alpha;
                }
            }
        }
    }
}

} } // namespace ulmblas, hpc

#endif // HPC_ULMBLAS_GESCAL_H
