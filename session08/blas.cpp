#include "blas.hpp"

//==============================================================================
// ulmBLAS
//==============================================================================

namespace ulmBLAS {

//------------------------------------------------------------------------------
// A -> B
//------------------------------------------------------------------------------

void
gecopy(const GeMatrix& A, GeMatrix& B)
{
    assert(A.m == B.m && A.n == B.n);
    for (std::size_t j=0; j<B.n; ++j) {
        for (std::size_t i=0; i<B.m; ++i) {
            B(i,j) = A(i,j);
        }
    }
}

//------------------------------------------------------------------------------
// Y <- Y + alpha*X
//------------------------------------------------------------------------------

void
geaxpy(double alpha, const GeMatrix& X, GeMatrix& Y)
{
    assert(X.m == Y.m && X.n == Y.n);
    if(alpha != 0.0) {
        for (std::size_t i=0; i<X.m; ++i) {
            for (std::size_t j=0; j<X.n; ++j) {
                Y(i,j) += alpha*X(i,j);
            }
        }
    }
}

//------------------------------------------------------------------------------
// A <- alpha * A
// if alpha = 0 initialize A with 0
//------------------------------------------------------------------------------

void
gescal(double alpha, GeMatrix& X)
{
    if (alpha!=1.0) {
        if(alpha != 0.0) {
            for (std::size_t i=0; i<X.m; ++i) {
                for (std::size_t j=0; j<X.n; ++j) {
                    X(i,j) *= alpha;
                }
            }
        } else {
            for (std::size_t i=0; i<X.m; ++i) {
                for (std::size_t j=0; j<X.n; ++j) {
                    X(i,j) = 0;
                }
            }
        }
    }
}

}   // namespace ulmBlas



//==============================================================================
// refColMajor
//==============================================================================

namespace refColMajor {

void
gemm(double alpha, const GeMatrix& A, const GeMatrix& B,
     double beta, GeMatrix& C)
{
    assert(C.m==A.m && C.n==B.n && A.n==B.m);
    ulmBLAS::gescal(beta, C);
    for (std::size_t j=0; j<C.n; ++j) {
        for (std::size_t l=0; l<A.n; ++l) {
            for (std::size_t i=0; i<C.m; ++i) {
                C(i,j) += alpha*A(i,l) * B(l,j);
            }
        }
    }
}

} // namespace refColMajor

/*
//==============================================================================
// simpleBuffer
//==============================================================================

#ifndef SIMPLE_PUFFER_M_C
#define SIMPLE_PUFFER_M_C 256
#endif

#ifndef SIMPLE_PUFFER_K_C
#define SIMPLE_PUFFER_K_C 256
#endif

#ifndef SIMPLE_PUFFER_N_C
#define SIMPLE_PUFFER_N_C 1024
#endif

namespace simpleBuffer {

double A_[SIMPLE_PUFFER_M_C*SIMPLE_PUFFER_K_C];
double B_[SIMPLE_PUFFER_K_C*SIMPLE_PUFFER_N_C];
double C_[SIMPLE_PUFFER_M_C*SIMPLE_PUFFER_N_C];

void
gemm(double alpha, const GeMatrix& A, const GeMatrix& B,
     double beta, GeMatrix& C)
{
    std::size_t M_C = SIMPLE_PUFFER_M_C;
    std::size_t N_C = SIMPLE_PUFFER_N_C;
    std::size_t K_C = SIMPLE_PUFFER_K_C;

    std::size_t mb = m / M_C;
    std::size_t nb = n / N_C;
    std::size_t kb = k / K_C;

    std::size_t mr = m % M_C;
    std::size_t nr = n % N_C;
    std::size_t kr = k % K_C;

    using ulmBLAS::geaxpy;
    using ulmBLAS::gecopy;
    using ulmBLAS::gescal;
    using refColMajor::gemm;

    gescal(m, n, beta, C, incRowC, incColC);

    for (std::size_t j=0; j<=nb; ++j) {
        std::size_t N = (j<nb || nr==0) ? N_C : nr;

        for (std::size_t l=0; l<=kb; ++l) {
            std::size_t K = (l<kb || kr==0) ? K_C : kr;

            gecopy(K, N,
                   &B[l*K_C*incRowB+j*N_C*incColB], incRowB, incColB,
                   B_, 1, K_C);

            for (std::size_t i=0; i<=mb; ++i) {
                std::size_t M = (i<mb || mr==0) ? M_C : mr;

                gecopy(M, K,
                       &A[i*M_C*incRowA+l*K_C*incColA], incRowA, incColA,
                       A_, 1, M_C);

                gemm(M, N, K,
                     1.0,
                     A_, 1, M_C,
                     B_, 1, K_C,
                     0.0,
                     C_, 1, M_C);

                geaxpy(M, N, alpha,
                       C_, 1, M_C,
                       &C[i*M_C*incRowC+j*N_C*incColC], incRowC, incColC);
            }
        }
    }
}

} // namespace simpleBuffer
*/
