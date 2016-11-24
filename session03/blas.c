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
    int i,j,l;
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
    int i,j,l;
    for(i=0; i<n; ++i) {
        for(j=0; j<m; ++j) {
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
    int i,j,l;
    for(i=0; i<m; ++i) {
        for(j=0; j<k; ++j) {
            double temp_A = A[i*incRowA + j*incColA];
            for(l=0; l<n; ++l) {
                double beta_ = (k==0) ? beta : 1.0;
                C[i*incRowC+l*incColC] *= beta_;
                C[i*incRowC+l*incColC] += alpha * (temp_A
                                                    * B[j*incRowB + l*incColB]);
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
    int i,j,l;
    for(i=0; i<n; ++i) {
        for(j=0; j<k; ++j) {
            double temp_B = B[j*incRowB + i*incColB];
            for(l=0; l<m; ++l) {
                double beta_ = (k==0) ? beta : 1.0;
                C[i*incRowC+l*incColC] *= beta_;
                C[i*incRowC+l*incColC] += alpha * ( A[l*incRowA+j*incColA]
                                                    * temp_B);
            }
        }
    }
}
