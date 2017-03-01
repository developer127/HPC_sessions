#include <stddef.h>

/** Matrix produkt C <- b*C + a*AB
    C m x n Matrix
    a,b scalars
    A m x k Matrix
    B k x n Matrix

    varinat 1   C row major     going through C once
                A row major
                B col major
*/
void gemm_variant1(size_t m, size_t n, size_t k,
                   double alpha,
                   const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                   const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                   double beta,
                   double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    size_t i,j,l;
    for(i=0; i<m; ++i) {
        for(j=0; j<n; ++j) {
            double temp = 0;
            for(l=0; l<k; ++l) {
                temp += A[i*incRowA+l*incColA] * B[l*incRowB + j*incColB];
            }
            C[i*incRowC+j*incColC] *= beta;
            C[i*incRowC+j*incColC] += temp * alpha;
        }
    }
}

/** varinat 2   C col major     goint trough C once
                A row major
                B col major
*/

void gemm_variant2(size_t m, size_t n, size_t k,
                   double alpha,
                   const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                   const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                   double beta,
                   double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    size_t i,j,l;
    for(j=0; j<n; ++j) {
        for(i=0; i<m; ++i) {
            double temp = 0;
            for(l=0; l<k; ++l) {
                temp += A[i*incRowA+l*incColA] * B[l*incRowB + j*incColB];
            }
            C[i*incRowC+j*incColC] *= beta;
            C[i*incRowC+j*incColC] += temp * alpha;
        }
    }
}

/** varinat 3   C row major
                A row major     going through A once
                B row major
*/

void gemm_variant3(size_t m, size_t n, size_t k,
                   double alpha,
                   const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                   const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                   double beta,
                   double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    size_t i,j,l;
    for(i=0; i<m; ++i) {
        for(l=0; l<k; ++l) {
            double beta_ = (l==0) ? beta : 1.0;
            double temp_A = A[i*incRowA+l*incColA];
            for(j=0; j<n; ++j) {
                C[i*incRowC+j*incColC] *= beta_;
                C[i*incRowC+j*incColC] += alpha * temp_A
                                        * B[l*incRowB+j*incColB];
            }
        }
    }
}

/** varinat 4   C col major
                A col major
                B col major     going through B once
*/

void gemm_variant4(size_t m, size_t n, size_t k,
                   double alpha,
                   const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                   const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                   double beta,
                   double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    size_t i,j,l;
    for(j=0; j<n; ++j) {
        for(l=0; l<k; ++l) {
            double beta_ = (l==0) ? beta : 1.0;
            double temp_B = B[l*incRowB+j*incColB];
            for(i=0; i<m; ++i) {
                C[i*incRowC+j*incColC] *= beta_;
                C[i*incRowC+j*incColC] += alpha * A[i*incRowA+l*incColA]
                                        * temp_B;
            }
        }
    }
}
