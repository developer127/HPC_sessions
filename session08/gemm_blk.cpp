#include <cmath>
#include <cstdio>
#include <sys/times.h>
#include <unistd.h>
#include <fmt/printf.hpp>

#ifndef COLMAJOR
#define COLMAJOR 1
#endif

#ifndef MAXDIM_M
#define MAXDIM_M    4000
#endif

#ifndef MAXDIM_N
#define MAXDIM_N    4000
#endif

#ifndef MAXDIM_K
#define MAXDIM_K    4000
#endif

#ifndef MIN_M
#define MIN_M   100
#endif

#ifndef MIN_N
#define MIN_N   100
#endif

#ifndef MIN_K
#define MIN_K   100
#endif

#ifndef MAX_M
#define MAX_M   7000
#endif

#ifndef MAX_N
#define MAX_N   7000
#endif

#ifndef MAX_K
#define MAX_K   7000
#endif

#ifndef INC_M
#define INC_M   100
#endif

#ifndef INC_N
#define INC_N   100
#endif

#ifndef INC_K
#define INC_K   100
#endif

#ifndef ALPHA
#define ALPHA   1.5
#endif

#ifndef BETA
#define BETA    1.5
#endif


//==============================================================================
// bench
//==============================================================================

namespace bench {

//------------------------------------------------------------------------------
// Auxiliary data for benchmarking
//------------------------------------------------------------------------------

double A[MAXDIM_M*MAXDIM_K];
double B[MAXDIM_K*MAXDIM_N];

double C1[MAXDIM_M*MAXDIM_N];
double C2[MAXDIM_M*MAXDIM_N];
double C3[MAXDIM_M*MAXDIM_N];

//------------------------------------------------------------------------------
// Auxiliary functions for benchmarking
//------------------------------------------------------------------------------

double
walltime()
{
   struct tms    ts;
   static double ClockTick=0.0;

   if (ClockTick==0.0) {
        ClockTick = 1.0 / ((double) sysconf(_SC_CLK_TCK));
   }
   return ((double) times(&ts)) * ClockTick;
}

void
initMatrix(std::size_t m, std::size_t n,
           double *A, std::size_t incRowA, std::size_t incColA)
{
    for (std::size_t j=0; j<n; ++j) {
        for (std::size_t i=0; i<m; ++i) {
            A[i*incRowA+j*incColA] = ((double)rand() - RAND_MAX/2)*200/RAND_MAX;
        }
    }
}

double
asumDiffMatrix(std::size_t m, std::size_t n,
               const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
               double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB)
{

    double asum = 0;

    for (std::size_t j=0; j<n; ++j) {
        for (std::size_t i=0; i<m; ++i) {
            asum += fabs(B[i*incRowB+j*incColB] - A[i*incRowA+j*incColA]);
        }
    }
    return asum;
}

void
printMatrixInMemory(std::size_t m, std::size_t n, const double *A)
{
    for(size_t i= 0; i< m*n; ++i) {
        fmt::printf("%5.1lf ", A[i]);
    }
    fmt::printf("\n");
}

void print_matrix(std::size_t m, std::size_t n, const double *A,
                  std::ptrdiff_t incRowA, std::ptrdiff_t incColA)
{
    for (std::size_t i = 0; i < m; ++i) {
        fmt::printf("  ");
        for (std::size_t j = 0; j < n; ++j) {
            fmt::printf(" %5.1lf", A[i*incRowA + j*incColA]);
        }
        fmt::printf("\n");
    }
}

} // namespace bench


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
gepack_A(std::size_t mc, std::size_t kc,
         const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
         double *p)
{
    std::size_t MR = BLOCKED_MR;
    std::size_t mp = (mc+MR-1) / MR;    //number of panels

    for(std::size_t j=0; j<kc; ++j) {
        for(std::size_t i=0; i<mp*MR; ++i) {
            std::size_t index = i/MR * kc * MR + j*MR + i%MR;
            p[index] = (i<mc)? A[i*incRowA+j*incColA] : 0;
        }
    }
}

void
gepack_B(std::size_t kc, std::size_t nc,
       const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
       double *p)
{
    std::size_t NR = BLOCKED_NR;
    std::size_t np = (nc+NR-1) / NR;    //number of panels

    for(std::size_t j=0; j<np*NR; ++j) {
        for(std::size_t i=0; i<kc; ++i) {
            std::size_t index = j/NR * kc * NR + i*NR + j%NR;
            p[index] = (j<nc)? B[i*incRowB+j*incColB] : 0;
        }
    }
}

void
ugemm(std::size_t kc, double alpha,
      const double *A, const double *B,
      double beta,
      double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC)
{
    double P[BLOCKED_MR*BLOCKED_NR];
    std::size_t MR = BLOCKED_MR;
    std::size_t NR = BLOCKED_NR;

    // ... FIXME
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

    // ... FIXME
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

int
main()
{
    using namespace bench;
    using namespace ulmBLAS;
    using namespace std;

    double packedA[20];
    double packedB[20];
    initMatrix(3, 5, A, 5, 1);
    print_matrix(3, 5, A, 5, 1);
    blocked::gepack_A(3, 5, A, 5, 1, packedA);
    printMatrixInMemory(4, 5, packedA);
    blocked::gepack_B(3, 5, A, 5, 1, packedB);
    printMatrixInMemory(3, 6, packedB);


    /*
    initMatrix(MAXDIM_M, MAXDIM_K, A, 1, MAXDIM_M);
    initMatrix(MAXDIM_K, MAXDIM_N, B, 1, MAXDIM_K);
    initMatrix(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M);
    gecopy(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M, C2, 1, MAXDIM_M);
    gecopy(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M, C3, 1, MAXDIM_M);

    // Header-Zeile für die Ausgabe
    fmt::printf("%5s %5s %5s ", "m", "n", "k");
    fmt::printf("%5s %5s ", "IRA", "ICA");
    fmt::printf("%5s %5s ", "IRB", "ICB");
    fmt::printf("%5s %5s ", "IRC", "ICC");
    fmt::printf("%20s %9s", "refColMajor: t", "MFLOPS");
    fmt::printf("%20s %9s %9s", "refSimpleBuffer: t", "MFLOPS", "diff");
    fmt::printf("\n");

    for (std::size_t m = MIN_M, n = MIN_N, k = MIN_K;
         m <=MAX_M && n <= MAX_N && k <= MAX_K;
         m += INC_M, n += INC_N, k += INC_K)
    {
        double t0, t1, diff;

        std::ptrdiff_t incRowA = (COLMAJOR) ? 1 : k;
        std::ptrdiff_t incColA = (COLMAJOR) ? m : 1;

        std::ptrdiff_t incRowB = (COLMAJOR) ? 1 : n;
        std::ptrdiff_t incColB = (COLMAJOR) ? k : 1;

        std::ptrdiff_t incRowC = (COLMAJOR) ? 1 : n;
        std::ptrdiff_t incColC = (COLMAJOR) ? m : 1;

        fmt::printf("%5zd %5zd %5zd ", m, n, k);
        fmt::printf("%5td %5td ", incRowA, incColA);
        fmt::printf("%5td %5td ", incRowB, incColB);
        fmt::printf("%5td %5td ", incRowC, incColC);

        t0 = walltime();
        refColMajor::gemm(m, n, k, ALPHA,
                          A, incRowA, incColA,
                          B, incRowB, incColB,
                          BETA,
                          C1, incRowC, incColC);
        t1 = walltime() - t0;
        fmt::printf("%20.4lf %9.2lf", t1, 2.*m/1000*n/1000*k/t1);

        t0 = walltime();
        simpleBuffer::gemm(m, n, k, ALPHA,
                           A, incRowA, incColA,
                           B, incRowB, incColB,
                           BETA,
                           C2, incRowC, incColC);
        t1 = walltime() - t0;
        diff = asumDiffMatrix(m, n,
                              C1, incRowC, incColC,
                              C2, incRowC, incColC);
        fmt::printf("%20.4lf %9.2lf %9.1e", t1, 2.*m/1000*n/1000*k/t1, diff);

        t0 = walltime();
        blocked::gemm(m, n, k, ALPHA,
                      A, incRowA, incColA,
                      B, incRowB, incColB,
                      BETA,
                      C3, incRowC, incColC);
        t1 = walltime() - t0;
        diff = asumDiffMatrix(m, n,
                              C1, incRowC, incColC,
                              C3, incRowC, incColC);
        fmt::printf("%20.4lf %9.2lf %9.1e", t1, 2.*m/1000*n/1000*k/t1, diff);
        fmt::printf("\n");
    }*/
    return 0;
}
