#ifndef INC_GECOPY_HPP
#define INC_GECOPY_HPP 1

namespace hpc { namespace ulmblas {

template <typename TA, typename TB, typename Size, typename Index>
void
gecopy(Size m, Size n,
       const TA& A, Index incRowA, Index incColA,
       TB&       B, Index incRowB, Index incColB)
{
    if(incColA < incRowA) {
        for(Size i=0; i<m; ++i) {
            for(Size j=0; j<n; ++j) {
                B[i*incRowB + j*incColB] = A[i*incRowA + j*incColA];
            }
        }
    } else {
        for(Size j=0; j<n; ++j) {
            for(Size i=0; i<m; ++i) {
                B[i*incRowB + j*incColB] = A[i*incRowA + j*incColA];
            }
        }

    }
}

}}      // namespace hpc::ulmblas
#endif  // INC_GECOPY_HPP
