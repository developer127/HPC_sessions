#ifndef INC_GEAXPY_HPP
#define INC_GEAXPY_HPP 1

namespace hpc { namespace ulmblas {

//------------------------------------------------------------------------------
// Y <- Y + alpha*X
//------------------------------------------------------------------------------

template <typename T, typename TX, typename TY, typename Size, typename Index>
void
geaxpy(Size m, Size n,
       T alpha,
       const TX *X, Index incRowX, Index incColX,
       TY       *Y, Index incRowY, Index incColY)
{
    if(alpha != T(0)) {
        if(incColX < incRowX) {
            for (Size i=0; i<m; ++i) {
                for (Size j=0; j<n; ++j) {
                    Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
                }
            }
        } else {
            for (Size j=0; j<n; ++j) {
                for (Size i=0; i<m; ++i) {
                    Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
                }
            }
        }
    }
}

}}      // namespace hpc::ulmblas

#endif  // INC_GEAXPY_HPP

