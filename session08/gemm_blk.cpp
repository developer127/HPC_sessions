#include <cmath>
#include <cstdio>
#include <sys/times.h>
#include <unistd.h>
#include <printf.hpp>
#include "benchmark.hpp"
#include "test_func.hpp"
#include "matrix_func.hpp"
#include "blas.hpp"

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

//------------------------------------------------------------------------------
// Auxiliary data for benchmarking
//------------------------------------------------------------------------------

double A[MAXDIM_M*MAXDIM_K];
double B[MAXDIM_K*MAXDIM_N];

double C1[MAXDIM_M*MAXDIM_N];
double C2[MAXDIM_M*MAXDIM_N];
double C3[MAXDIM_M*MAXDIM_N];


int
main()
{
    using namespace bench;
    using namespace ulmBLAS;
    using namespace std;
    using namespace test;

    double A[3*3];
    double B[3*3];
    double C_1[3*3];
    double C_2[3*3];

    initMatrix(3, 3, A, 1, 3);
    initMatrix(3, 4, B, 3, 1);
    initMatrix(3, 3, C_1, 3, 1);
    gecopy(3, 3,
           C_1, 3, 1,
           C_2, 3, 1);

    refColMajor::gemm(3, 3, 3, 1,
                      A, 1, 3,
                      B, 3, 1,
                      1,
                      C1, 3, 1);

    double A_[4*3];
    double B_[3*4];
    blocked::pack_A(3,3,
                    A, 1, 3,
                    A_);

    printMatrix(3,3, A, 1, 3);
    test::printGeMatrixInMemory(3, 3, A);
    test::printGeMatrixInMemory(4, 3, A_);

    blocked::pack_B(3,3,
                    B, 3, 1,
                    B_);

    printMatrix(3,3, B, 3, 1);
    test::printGeMatrixInMemory(3, 3, B);
    test::printGeMatrixInMemory(3, 4, B_);

    blocked::mgemm(3, 3, 1, 1,
                   A_, B_,
                   1,
                   C2, 3, 1);

    double diff = asumDiffMatrix(3, 3,
                                 C_1, 3, 1,
                                 C_2, 3, 1);

    fmt::printf("Differenze: %lf\n", diff);
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
    }
    return 0;*/
}
