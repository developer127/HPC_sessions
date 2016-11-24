#include <cstdio>
#include "bench.hpp"
#include "ulmblas.hpp"
#include "gemm_refcolmajor.hpp"
#include "gemm_blocked.hpp"

#ifndef COLMAJOR
#define COLMAJOR 1
#endif

#ifndef MAXDIM_M
#define MAXDIM_M    7000
#endif

#ifndef MAXDIM_N
#define MAXDIM_N    7000
#endif

#ifndef MAXDIM_K
#define MAXDIM_K    7000
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

int
main()
{

    typedef double      T;

    T *A  = new T[MAXDIM_M*MAXDIM_K];
    T *B  = new T[MAXDIM_K*MAXDIM_N];
    T *C1 = new T[MAXDIM_M*MAXDIM_N];
    T *C2 = new T[MAXDIM_M*MAXDIM_N];

    bench::initGeMatrix(MAXDIM_M, MAXDIM_K, A, 1, MAXDIM_M);
    bench::initGeMatrix(MAXDIM_K, MAXDIM_N, B, 1, MAXDIM_K);
    bench::initGeMatrix(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M);
    ulmBLAS::gecopy(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M, C2, 1, MAXDIM_M);

    // Header-Zeile fuer die Ausgabe
    printf("%5s %5s %5s ", "m", "n", "k");
    printf("%5s %5s ", "IRA", "ICA");
    printf("%5s %5s ", "IRB", "ICB");
    printf("%5s %5s ", "IRC", "ICC");
    printf("%20s %9s", "refColMajor: t", "MFLOPS");
    printf("%20s %9s %9s", "blocked GEMM: t", "MFLOPS", "diff");
    printf("\n");

    bench::WallTime<double> wallTime;

    for (long m = MIN_M, n = MIN_N, k = MIN_K;
         m <=MAX_M && n <= MAX_N && k <= MAX_K;
         m += INC_M, n += INC_N, k += INC_K)
    {
        double t, diff;

        long incRowA = (COLMAJOR) ? 1 : k;
        long incColA = (COLMAJOR) ? m : 1;

        long incRowB = (COLMAJOR) ? 1 : n;
        long incColB = (COLMAJOR) ? k : 1;

        long incRowC = (COLMAJOR) ? 1 : n;
        long incColC = (COLMAJOR) ? m : 1;

        printf("%5ld %5ld %5ld ", m, n, k);
        printf("%5ld %5ld ", incRowA, incColA);
        printf("%5ld %5ld ", incRowB, incColB);
        printf("%5ld %5ld ", incRowC, incColC);

        wallTime.tic();
        refColMajor::gemm(m, n, k, ALPHA,
                          A, incRowA, incColA,
                          B, incRowB, incColB,
                          BETA,
                          C1, incRowC, incColC);
        t = wallTime.toc();
        printf("%20.4lf %9.2lf", t, 2.*m/1000*n/1000*k/t);

        wallTime.tic();
        blocked::gemm(m, n, k, ALPHA,
                      A, incRowA, incColA,
                      B, incRowB, incColB,
                      BETA,
                      C2, incRowC, incColC);
        t = wallTime.toc();
        diff = bench::asumDiffGeMatrix(m, n,
                                       C1, incRowC, incColC,
                                       C2, incRowC, incColC)/(m*n);
        printf("%20.4lf %9.2lf %9.1e", t, 2.*m/1000*n/1000*k/t, diff);
        printf("\n");
    }

    delete [] A;
    delete [] B;
    delete [] C1;
    delete [] C2;
}
