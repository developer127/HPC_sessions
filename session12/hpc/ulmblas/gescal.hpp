#ifndef INC_GESCAL_HPP
#define INC_GESCAL_HPP 1

namespace hpc { namespace ulmblas {


template <typename T, typename TX, typename Size, typename Index>
void
gescal(Size m, Size n, T alpha, TX X, Index incRowX, Index incColX)
{
    if (alpha != T(1)) {
        if (alpha != T(0)) {
            if(incColX < incRowX) {
                for (Size i=0; i<m; ++i) {
                    for (Size j=0; j<n; ++j) {
                        X[i*incRowX+j*incColX] *= alpha;
                    }
                }
            } else {
                for (Size j=0; j<n; ++j) {
                    for (Size i=0; i<m; ++i) {
                        X[i*incRowX+j*incColX] *= alpha;
                    }
                }
            }
        } else {
            if(incColX < incRowX) {
                for (Size i=0; i<m; ++i) {
                    for (Size j=0; j<n; ++j) {
                        X[i*incRowX+j*incColX] = T(0);
                    }
                }
            } else {
                for (Size j=0; j<n; ++j) {
                    for (Size i=0; i<m; ++i) {
                        X[i*incRowX+j*incColX] = T(0);
                    }
                }

            }
        }
    }

}

}}      // namespace hpc::ulmblas

#endif  // INC_GESCAL_HPP
