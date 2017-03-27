#include "blas.hpp"

//==============================================================================
// ulmBLAS
//==============================================================================

namespace ulmBLAS {

//------------------------------------------------------------------------------
// A <- B
//------------------------------------------------------------------------------

void
gecopy(std::size_t m, std::size_t n,
       const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
       double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB)
{
    for (std::size_t j=0; j<n; ++j) {
        for (std::size_t i=0; i<m; ++i) {
            B[i*incRowB+j*incColB] = A[i*incRowA+j*incColA];
        }
    }
}

//------------------------------------------------------------------------------
// Y <- Y + alpha*Yä#
//------------------------------------------------------------------------------

void
geaxpy(std::size_t m, std::size_t n,
       double alpha,
       const double *X, std::ptrdiff_t incRowX, std::ptrdiff_t incColX,
       double       *Y, std::ptrdiff_t incRowY, std::ptrdiff_t incColY)
{
    for (std::size_t i=0; i<m; ++i) {
        for (std::size_t j=0; j<n; ++j) {
            Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
        }
    }
}

//------------------------------------------------------------------------------
// A <- alpha * A
//------------------------------------------------------------------------------

void
gescal(std::size_t m, std::size_t n,
       double alpha,
       double *X, std::ptrdiff_t incRowX, std::ptrdiff_t incColX)
{
    if (alpha!=1.0) {
        for (std::size_t i=0; i<m; ++i) {
            for (std::size_t j=0; j<n; ++j) {
                X[i*incRowX+j*incColX] *= alpha;
            }
        }
    }
}

//------------------------------------------------------------------------------
// C <- alpha x * y     rank1 update
//------------------------------------------------------------------------------

void
ger(std::size_t m, std::size_t n,
    double alpha,
    const double *x, std::ptrdiff_t incX,
    const double *y, std::ptrdiff_t incY,
    double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC)
{
    if (alpha!=1.0) {
        if (incRowC < incColC) {    //colMajor
            for (std::size_t j=0; j<n; ++j) {
                double tmp_y = y[j*incY];
                for (std::size_t i=0; i<m; ++i) {
                    C[i*incRowC+j*incColC] += alpha * x[i*incX] * tmp_y;
                }
            }
        } else {                    //rowMajor
            for (std::size_t i=0; i<m; ++i) {
                double tmp_x = x[i*incX];
                for (std::size_t j=0; j<n; ++j) {
                    C[i*incRowC+j*incColC] += alpha * tmp_x * y[j*incY];
                }
            }
        }
    }
}


} // namespace ulmBLAS

//==============================================================================
// refColMajor
//==============================================================================

namespace refColMajor {

void
gemm(std::size_t m, std::size_t n, std::size_t k,
     double alpha,
     const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
     const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
     double beta,
     double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC)
{
    ulmBLAS::gescal(m, n, beta, C, incRowC, incColC);
    for (std::size_t j=0; j<n; ++j) {
        for (std::size_t l=0; l<k; ++l) {
            for (std::size_t i=0; i<m; ++i) {
                C[i*incRowC+j*incColC] += alpha*A[i*incRowA+l*incColA]
                                               *B[l*incRowB+j*incColB];
            }
        }
    }
}

} // namespace refColMajor



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
gemm(std::size_t m, std::size_t n, std::size_t k,
     double alpha,
     const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
     const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
     double beta,
     double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC)
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


//==============================================================================
// blocked
//==============================================================================


#ifndef BLOCKED_MC
#define BLOCKED_MC 256
#endif

#ifndef BLOCKED_KC
#define BLOCKED_KC 256
#endif

#ifndef BLOCKED_NC
#define BLOCKED_NC 1024
#endif

#ifndef BLOCKED_MR
#define BLOCKED_MR 2
#endif

#ifndef BLOCKED_NR
#define BLOCKED_NR 8
#endif

namespace blocked {

double A_[BLOCKED_MC*BLOCKED_KC];
double B_[BLOCKED_KC*BLOCKED_NC];

void
pack_A(std::size_t mc, std::size_t kc,
       const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
       double *p)
{
    std::size_t MR = BLOCKED_MR;
    std::size_t mp = (mc+MR-1) / MR;        // nof rowPanels in A

    if (incRowA < incColA) {                // colMajor
        for (std::size_t j=0; j<kc; ++j) {
            for (std::size_t i=0; i< mp*MR; ++i) {
                std::size_t l = i/MR * kc * MR + j*MR + i%MR;
                p[l] = (i<mc) ? A[i*incRowA+j*incColA] : 0;
            }
        }
    } else {
        for (std::size_t i=0; i< mp*MR; ++i) {
            for (std::size_t j=0; j<kc; ++j) {
                std::size_t l = i/MR * kc * MR + j*MR + i%MR;
                p[l] = (i<mc) ? A[i*incRowA+j*incColA] : 0;
            }
        }
    }
}

void
pack_B(std::size_t kc, std::size_t nc,
       const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
       double *p)
{
    std::size_t NR = BLOCKED_NR;
    std::size_t np = (nc+NR-1) / NR;        // nof colPanels in B

    if (incRowB < incColB) {                // colMajor
        for (std::size_t j=0; j<np*NR; ++j) {
            for (std::size_t i=0; i<kc ; ++i) {
                std::size_t l = j/NR * kc * NR + i*NR + j%NR;
                p[l] = (j<nc) ? B[i*incRowB+j*incColB] : 0;
            }
        }
    } else {
        for (std::size_t i=0; i<kc ; ++i) {
            for (std::size_t j=0; j<np*NR; ++j) {
                std::size_t l = j/NR * kc * NR + i*NR + j%NR;
                p[l] = (j<nc) ? B[i*incRowB+j*incColB] : 0;
            }
        }
    }
}

void
ugemm(std::size_t kc, double alpha,
      const double *A, const double *B,
      double beta,
      double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC)
{
    for (std::size_t k=0; k<kc; ++k) {
        ulmBLAS::ger(BLOCKED_MR, BLOCKED_NR, alpha,
                     &A[k*BLOCKED_MR], 1,
                     &B[k*BLOCKED_NR], 1,
                     C, incRowC, incColC);
    }
}

void
mgemm(std::size_t mc, std::size_t nc, std::size_t kc,
      double alpha,
      const double *A, const double *B,
      double beta,
      double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC)
{
    double C_[BLOCKED_MR*BLOCKED_NR];
    std::size_t MR = BLOCKED_MR;
    std::size_t NR = BLOCKED_NR;

    std::size_t mp  = (mc+MR-1) / MR;
    std::size_t np  = (nc+NR-1) / NR;
    std::size_t mr_ = mc % MR;
    std::size_t nr_ = nc % NR;

    for (std::size_t j=0; j< np; ++j) {
        std::size_t nr = (j<np-1 || nr_==0) ? NR : nr_;
        for (std::size_t i=0; i< mp; ++i) {
            std::size_t mr = (i<mp-1 || mr_==0) ? MR : mr_;
            if (mr==MR && nr==NR) {
                ugemm(kc, alpha,
                      &A[i*kc*MR], &B[j*kc*NR], beta,
                      &C[i*MR*incRowC + j*NR*incColC], incRowC, incColC);
            } else {
                ulmBLAS::gecopy(mr, nr,
                                &C[i*MR*incRowC + j*NR*incColC],
                                incRowC, incColC,
                                C_, 1, MR);      // C_ is colMajor

                ugemm(kc, alpha,
                      &A[i*kc*MR], &B[j*kc*NR], beta,
                      C_, incRowC, incColC);

                ulmBLAS::gecopy(mr, nr,
                                C_, 1, MR,
                                &C[i*MR*incRowC + j*NR*incColC],
                                incRowC, incColC);
            }
        }
    }
}

void
gemm(std::size_t m, std::size_t n, std::size_t k,
     double alpha,
     const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
     const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
     double beta,
     double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC)
{
    std::size_t MC = BLOCKED_MC;
    std::size_t NC = BLOCKED_NC;
    std::size_t KC = BLOCKED_KC;

    std::size_t mb = (m+MC-1) / MC;
    std::size_t nb = (n+NC-1) / NC;
    std::size_t kb = (k+KC-1) / KC;

    std::size_t mc_ = m % MC;
    std::size_t nc_ = n % NC;
    std::size_t kc_ = k % KC;

    // ... FIXME
}

} // namespace blocked

