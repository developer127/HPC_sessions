#ifndef INC_GEMM_REFCOLMAJOR_H
#define INC_GEMM_REFCOLMAJOR_H 1

namespace refColMajor {

template <typename T, typename Size, typename Index>
void
gemm(Size m, Size n, Index k,
    T alpha,
    const T *A, Index incRowA, Index incColA,
    const T *B, Index incRowB, Index incColB,
    T beta,
    T *C, Index incRowC, Index incColC)
{
    ulmBLAS::gescal(m, n, beta, C, incRowC, incColC);
    if(alpha == T(0)) {
        return;
    }
    for (Size j=0; j<n; ++j) {
        for (Size l=0; l<k; ++l) {
            for (Size i=0; i<m; ++i) {
                C[i*incRowC+j*incColC] += alpha*A[i*incRowA+l*incColA]
                                               *B[l*incRowB+j*incColB];
            }
        }
    }
}

} // namespace refColMajor

#endif // INC_GEMM_REFCOLMAJOR_H
