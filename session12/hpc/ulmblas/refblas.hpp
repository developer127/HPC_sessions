#ifndef INC_REFBLAS_HPP
#define INC_REFBLAS_HPP 1


namespace refblas {

template <typename Alpha, typename TA, typename TB, typename Beta, typename TC,
          typename Size, typename Index>
void
gemm(Size m, Size n, Size k,
     Alpha alpha,
     const TA *A, Index incRowA, Index incColA,
     const TB *B, Index incRowB, Index incColB,
     Beta beta,
     TC *C, Index incRowC, Index incColC)
{
    hpc::ulmblas::gescal(m, n, beta, C, incRowC, incColC);
    for (Index j=0; j<n; ++j) {
        for (Index l=0; l<k; ++l) {
            for (Index i=0; i<m; ++i) {
                C[i*incRowC+j*incColC] += alpha*A[i*incRowA+l*incColA]
                                               *B[l*incRowB+j*incColB];
            }
        }
    }
}

}       // namespace refblas

#endif  // INC_REFBLAS_HPP
