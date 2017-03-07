#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
#include <stddef.h>

#ifndef COLMAJOR
#define COLMAJOR 0
#else
#undef  COLMAJOR
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
#define MIN_M   200
#endif

#ifndef MIN_N
#define MIN_N   200
#endif

#ifndef MIN_K
#define MIN_K   200
#endif

#ifndef MAX_M
#define MAX_M   4000
#endif

#ifndef MAX_N
#define MAX_N   4000
#endif

#ifndef MAX_K
#define MAX_K   4000
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
#define ALPHA   1.0
#endif

#ifndef BETA
#define BETA    1.5
#endif

#ifndef M_C
#define M_C     100
#endif

#ifndef K_C
#define K_C     100
#endif

#ifndef N_C
#define N_C     100
#endif



double A[MAXDIM_M*MAXDIM_K];
double B[MAXDIM_K*MAXDIM_N];

double C1[MAXDIM_M*MAXDIM_N];
double C2[MAXDIM_M*MAXDIM_N];
double C3[MAXDIM_M*MAXDIM_N];
double C4[MAXDIM_M*MAXDIM_N];

//Buffer Matrix
double A_[M_C*K_C];
double B_[K_C*N_C];
double C_[M_C*N_C];

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
initMatrix(size_t m, size_t n, double *A, ptrdiff_t incRowA, ptrdiff_t incColA)
{
    for (size_t j=0; j<n; ++j) {
        for (size_t i=0; i<m; ++i) {
            A[i*incRowA+j*incColA] = ((double)rand() - RAND_MAX/2)*200/RAND_MAX;
        }
    }
}

// faster if C is col_major
void
initZero(size_t m, size_t n, double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    for(size_t j = 0; j<n; ++j) {
        for(size_t i=0; i<m; ++i) {
            C[i*incRowC+j*incColC] = 0;
        }
    }
}

void
dgecopy(size_t m, size_t n,
        const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
        double *B, ptrdiff_t incRowB, ptrdiff_t incColB)
{
    for (size_t j=0; j<n; ++j) {
        for (size_t i=0; i<m; ++i) {
            B[i*incRowB+j*incColB] = A[i*incRowA+j*incColA];
        }
    }
}

void
dgeaxpy(size_t m, size_t n,
        double alpha,
        const double *X, ptrdiff_t incRowX, ptrdiff_t incColX,
        double       *Y, ptrdiff_t incRowY, ptrdiff_t incColY)
{
    for (size_t i=0; i<m; ++i) {
        for (size_t j=0; j<n; ++j) {
            Y[i*incRowY+j*incColY] += alpha*X[i*incRowX+j*incColX];
        }
    }
}

double
asumMatrix(size_t m, size_t n,
           const double *A, ptrdiff_t incRowA, ptrdiff_t incColA)
{
    double asum = 0;

    for (size_t j=0; j<n; ++j) {
        for (size_t i=0; i<m; ++i) {
            asum += fabs(A[i*incRowA+j*incColA]);
        }
    }
    return asum;
}


double
asumDiffMatrix(size_t m, size_t n,
               const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
               const double *B, ptrdiff_t incRowB, ptrdiff_t incColB)
{
    size_t i, j;

    double asum = 0;

    for (j=0; j<n; ++j) {
        for (i=0; i<m; ++i) {
            asum += fabs(B[i*incRowB+j*incColB] - A[i*incRowA+j*incColA]);
        }
    }
    return asum;
}

#define MAX2(x,y)     ((x)>(y)) ? (x) : (y)

double
estimateGemmResidual(size_t m, size_t n, size_t k,
                     double alpha,
                     const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
                     const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
                     double beta,
                     const double *C0, ptrdiff_t incRowC0, ptrdiff_t incColC0,
                     const double *C1, ptrdiff_t incRowC1, ptrdiff_t incColC1)
{
    double aNorm = fabs(alpha)*asumMatrix(m, k, A, incRowA, incColA);
    double bNorm = asumMatrix(k, n, B, incRowB, incColB);
    double cNorm = fabs(beta)*asumMatrix(k, n, B, incRowB, incColB);
    double diff  = asumDiffMatrix(m, n,
                                  C0, incRowC0, incColC0,
                                  C1, incRowC1, incColC1);
    return diff/(aNorm*bNorm*cNorm*DBL_EPSILON*MAX2(MAX2(m, n), k));
}



//------------------------------------------------------------------------------
// GEMM simple block    schoolish version
//------------------------------------------------------------------------------
void
gemm_block(const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
           const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
           double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    for (size_t i=0; i<M_C; ++i) {
        for (size_t j=0; j<N_C; ++j) {
            double dot = 0;
            for(size_t l=0; l<K_C; ++l) {
                dot += A[i*incRowA+l*incColA]
                    * B[l*incRowB+j*incColB];
            }
            C[i*incRowC+j*incColC] += dot;
        }
    }
}

//------------------------------------------------------------------------------
// 1. GEMM Variant: School method
//------------------------------------------------------------------------------
void
dgemm_var1(size_t m, size_t n, size_t k, double alpha,
           const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
           const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
           double beta,
           double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    // Compute:  C <- beta*C
    if (beta!=1.0) {
        if (beta!=0) {
            for (size_t i=0; i<m; ++i) {
                for (size_t j=0; j<n; ++j) {
                    C[i*incRowC+j*incColC] *=beta;
               }
            }
        } else {
            for (size_t i=0; i<m; ++i) {
                for (size_t j=0; j<n; ++j) {
                    C[i*incRowC+j*incColC] = 0;
                }
            }
        }
    }
    // Compute: C <- C + alpha*A*B
    if(alpha != 0) {
        for (size_t i=0; i<m; ++i) {
            for (size_t j=0; j<n; ++j) {
                double dot = 0;
                for(size_t l=0; l<k; ++l) {
                    dot += A[i*incRowA+l*incColA]
                        * B[l*incRowB+j*incColB];
                }
                C[i*incRowC+j*incColC] += alpha * dot;
            }
        }
    }
}


//------------------------------------------------------------------------------
// 2. GEMM Variant: traverse row wise
//------------------------------------------------------------------------------
void
dgemm_var2(size_t m, size_t n, size_t k, double alpha,
           const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
           const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
           double beta,
           double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    // Compute:  C <- beta*C
    if (beta!=1.0) {
        if (beta!=0) {
            for (size_t i=0; i<m; ++i) {
                for (size_t j=0; j<n; ++j) {
                    C[i*incRowC+j*incColC] *=beta;
               }
            }
        } else {
            for (size_t i=0; i<m; ++i) {
                for (size_t j=0; j<n; ++j) {
                    C[i*incRowC+j*incColC] = 0;
                }
            }
        }
    }
    // Compute: C <- C + alpha*A*B
    if(alpha == 0) {
        return;
    }
    for(size_t i=0; i<m; ++i) {
        for(size_t l=0; l<k; ++l) {
            double temp_A = A[i*incRowA + l*incColA];
            for(size_t j=0; j<n; ++j) {
                C[i*incRowC+j*incColC] += alpha
                                        * (temp_A * B[l*incRowB + j*incColB]);
            }
        }
    }
}


//------------------------------------------------------------------------------
// 3. GEMM Variant: traverse column wise
//------------------------------------------------------------------------------
void
dgemm_var3(size_t m, size_t n, size_t k, double alpha,
           const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
           const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
           double beta,
           double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    // Compute:  C <- beta*C
    if (beta!=1.0) {
        if (beta!=0) {
            for (size_t i=0; i<m; ++i) {
                for (size_t j=0; j<n; ++j) {
                    C[i*incRowC+j*incColC] *=beta;
               }
            }
        } else {
            for (size_t i=0; i<m; ++i) {
                for (size_t j=0; j<n; ++j) {
                    C[i*incRowC+j*incColC] = 0;
                }
            }
        }
    }
    // Compute: C <- C + alpha*A*B
    if(alpha == 0) {
        return;
    }
    for(size_t j=0; j<n; ++j) {
        for(size_t l=0; l<k; ++l) {
            double temp_B = B[l*incRowB + j*incColB];
            for(size_t i=0; i<m; ++i) {
                C[i*incRowC+j*incColC] += alpha
                                        * A[i*incRowA+l*incColA]* temp_B;
            }
        }
    }
}


//------------------------------------------------------------------------------
// 4. GEMM Variant: traverse best
//------------------------------------------------------------------------------
void
dgemm_var4(size_t m, size_t n, size_t k, double alpha,
           const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
           const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
           double beta,
           double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    // Compute:  C <- beta*C
    if (beta!=1.0) {
        if (beta!=0) {
            if(incRowC>incColC){    // row major
                for (size_t i=0; i<m; ++i) {
                    for (size_t j=0; j<n; ++j) {
                        C[i*incRowC+j*incColC] *=beta;
                    }
                }
            } else {                // col major
                for (size_t j=0; j<m; ++j) {
                    for (size_t i=0; i<n; ++i) {
                        C[i*incRowC+j*incColC] *=beta;
                    }
                }

            }
        } else {
            if(incRowC>incColC){    // row major
                for (size_t i=0; i<m; ++i) {
                    for (size_t j=0; j<n; ++j) {
                        C[i*incRowC+j*incColC] = 0;
                    }
                }
            } else {                // col-major
                for (size_t j=0; j<m; ++j) {
                    for (size_t i=0; i<n; ++i) {
                        C[i*incRowC+j*incColC] = 0;
                    }
                }

            }
        }
    }
    // Compute: C <- C + alpha*A*B
    if(alpha == 0) {
        return;
    }
    if(incRowC>incColC){    // row major
        for(size_t i=0; i<m; ++i) {
            for(size_t l=0; l<k; ++l) {
                double temp_A = A[i*incRowA + l*incColA];
                for(size_t j=0; j<n; ++j) {
                    C[i*incRowC+j*incColC] += alpha
                                            * (temp_A * B[l*incRowB + j*incColB]);
                }
            }
        }
    } else {                // col-major
        for(size_t j=0; j<n; ++j) {
            for(size_t l=0; l<k; ++l) {
                double temp_B = B[l*incRowB + j*incColB];
                for(size_t i=0; i<m; ++i) {
                    C[i*incRowC+j*incColC] += alpha
                                            * A[i*incRowA+l*incColA]* temp_B;
                }
            }
        }
    }
}


//------------------------------------------------------------------------------
// 5. GEMM Variant: bockwise buffered
//------------------------------------------------------------------------------
void
dgemm_var5(size_t m, size_t n, size_t k, double alpha,
           const double *A, ptrdiff_t incRowA, ptrdiff_t incColA,
           const double *B, ptrdiff_t incRowB, ptrdiff_t incColB,
           double beta,
           double *C, ptrdiff_t incRowC, ptrdiff_t incColC)
{
    // Compute:  C <- beta*C
    if (beta!=1.0) {
        if (beta!=0) {
            if(incRowC>incColC){    // row major
                for (size_t i=0; i<m; ++i) {
                    for (size_t j=0; j<n; ++j) {
                        C[i*incRowC+j*incColC] *=beta;
                    }
                }
            } else {                // col major
                for (size_t j=0; j<m; ++j) {
                    for (size_t i=0; i<n; ++i) {
                        C[i*incRowC+j*incColC] *=beta;
                    }
                }

            }
        } else {
            if(incRowC>incColC){    // row major
                for (size_t i=0; i<m; ++i) {
                    for (size_t j=0; j<n; ++j) {
                        C[i*incRowC+j*incColC] = 0;
                    }
                }
            } else {                // col-major
                for (size_t j=0; j<m; ++j) {
                    for (size_t i=0; i<n; ++i) {
                        C[i*incRowC+j*incColC] = 0;
                    }
                }

            }
        }
    }
    // Compute: C <- C + alpha*A*B
    if(alpha == 0) {
        return;
    }
    // Anzahl der in die Matritzen passenden Blöcke
    size_t mb = m/M_C;
    size_t nb = n/N_C;
    size_t kb = k/K_C;

    // The remainder
    size_t mr = m % M_C;
    size_t nr = n % N_C;
    size_t kr = k % K_C;

    for(size_t i=0; i<=mb; ++i) {
        size_t M = (i<mb || mr==0)? M_C : mr;
        for(size_t j=0; j<=nb; ++j) {
            size_t N = (j<nb || nr==0)? N_C : nr;
            // Buffer C_ col_major
            initZero(M, N,
                     C_, 1, M_C);
            for(size_t l=0; l<=kb; ++l) {
                size_t K = (l<kb || kr==0)? K_C : kr;
                // Buffer A_ row_major
                dgecopy(M, K,
                        &A[i*incRowA*M_C+l*incColA*K_C], incRowA, incColA,
                        A_, K_C, 1);
                // Buffer B_ col_major
                dgecopy(K, N,
                        &B[l*incRowB*K_C+j*incColB*N_C], incRowB, incColB,
                        B_, 1, K_C);
                dgemm_var4(M, N, K, alpha,
                           A_, K_C, 1,
                           B_, 1, K_C,
                           1,
                           C_, 1, M_C);
            }
            dgeaxpy(M, N, 1,
                    C_, 1, M_C,
                    &C[i*incRowC*M_C+j*incColC*N_C], incRowC, incColC);
        }
    }
}


int
main()
{
    initMatrix(MAXDIM_M, MAXDIM_K, A, 1, MAXDIM_M);
    initMatrix(MAXDIM_K, MAXDIM_N, B, 1, MAXDIM_K);
    initMatrix(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M);
    dgecopy(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M, C2, 1, MAXDIM_M);
    dgecopy(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M, C3, 1, MAXDIM_M);
    dgecopy(MAXDIM_M, MAXDIM_N, C1, 1, MAXDIM_M, C4, 1, MAXDIM_M);

    // FIXME: Print a header for the output
    printf("%*c", 19, ' ');
    printf("%16s ", "GEMM Variant 1");
    //printf("%*c", 2, ' ');
    printf("%16s ", "GEMM Variant 2");
    printf("%*c", 10, ' ');
    printf("%16s ", "GEMM BEST");
    printf("%*c", 10, ' ');
    printf("%16s ", "GEMM BLOCKED");
    printf("\n");

    printf("%5s %5s %5s ", "m", "n", "k");
    printf("%7s %8s ", "t", "mflops");
    printf("%7s %8s %10s", "t", "mflops", "res");
    printf("%7s %8s %10s", "t", "mflops", "res");
    printf("%7s %8s %10s", "t", "mflops", "res");
    printf("\n");

    for (size_t m = MIN_M, n = MIN_N, k = MIN_K;
         m <=MAX_M && n <= MAX_N && k <= MAX_K;
         m += INC_M, n += INC_N, k += INC_K)
    {
        double t, mflops, res;

        size_t incRowA = (COLMAJOR) ? 1 : k;
        size_t incColA = (COLMAJOR) ? m : 1;

        size_t incRowB = (COLMAJOR) ? 1 : n;
        size_t incColB = (COLMAJOR) ? k : 1;

        size_t incRowC = (COLMAJOR) ? 1 : n;
        size_t incColC = (COLMAJOR) ? m : 1;

        printf("%5zu %5zu %5zu ", m, n, k);

        // Bench GEMM variant 1
        t = walltime();
        dgemm_var1(m, n, k, ALPHA,
                   A, incRowA, incColA,
                   B, incRowB, incColB,
                   BETA,
                   C1, incRowC, incColC);
        t = walltime() - t;
        mflops = 2.*m/1000*n/1000*k/t;

        printf("%7.2lf %8.2lf ", t, mflops);

        // Bench GEMM variant 2
        t = walltime();
        dgemm_var2(m, n, k, ALPHA,
                   A, incRowA, incColA,
                   B, incRowB, incColB,
                   BETA,
                   C2, incRowC, incColC);
        t = walltime() - t;
        mflops = 2.*m/1000*n/1000*k/t;
        res = estimateGemmResidual(m, n, k, ALPHA,
                                   A, incRowA, incColA,
                                   B, incRowB, incColB,
                                   BETA,
                                   C1, incRowC, incColC,
                                   C2, incRowC, incColC);

        printf("%7.2lf %8.2lf %10.2e", t, mflops, res);

        // Bench GEMM variant best choice
        t = walltime();
        dgemm_var4(m, n, k, ALPHA,
                   A, incRowA, incColA,
                   B, incRowB, incColB,
                   BETA,
                   C3, incRowC, incColC);
        t = walltime() - t;
        mflops = 2.*m/1000*n/1000*k/t;
        res = estimateGemmResidual(m, n, k, ALPHA,
                                   A, incRowA, incColA,
                                   B, incRowB, incColB,
                                   BETA,
                                   C1, incRowC, incColC,
                                   C3, incRowC, incColC);

        printf("%7.2lf %8.2lf %10.2e", t, mflops, res);

        // Bench GEMM variant BLOCKWISE
        t = walltime();
        dgemm_var5(m, n, k, ALPHA,
                   A, incRowA, incColA,
                   B, incRowB, incColB,
                   BETA,
                   C4, incRowC, incColC);
        t = walltime() - t;
        mflops = 2.*m/1000*n/1000*k/t;
        res = estimateGemmResidual(m, n, k, ALPHA,
                                   A, incRowA, incColA,
                                   B, incRowB, incColB,
                                   BETA,
                                   C1, incRowC, incColC,
                                   C4, incRowC, incColC);

        printf("%7.2lf %8.2lf %10.2e", t, mflops, res);

        printf("\n");
    }
}
