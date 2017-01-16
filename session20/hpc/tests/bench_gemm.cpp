#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <fmt/printf.hpp>
#include <hpc/aux/primitive_type.h>
#include <hpc/aux/walltime.h>
#include <hpc/mt/thread_pool.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/asum.h>
#include <hpc/matvec/axpy.h>
#include <hpc/matvec/copy.h>
#include <hpc/matvec/scal.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/isgematrix.h>
#include <hpc/matvec/mm.h>
#include <hpc/matvec/print.h>

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

template <typename MA>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value,
         void>::type
randomInit(MA &A)
{
    typedef typename MA::ElementType  T;
    typedef typename MA::Index        Index;

    std::random_device                  random;
    std::mt19937                        mt(random());
    std::uniform_real_distribution<T>   uniform(-100,100);

    auto func = [&](T &val, Index i, Index j) mutable -> void
                {
                    val = uniform(mt);
                };

    hpc::matvec::apply(A, func);
}

template <typename MA>
typename std::enable_if<hpc::matvec::IsComplexGeMatrix<MA>::value,
         void>::type
randomInit(MA &A)
{
    typedef typename MA::ElementType                    ET;
    typedef typename hpc::aux::PrimitiveType<ET>::Type  PT;
    typedef typename MA::Index                          Index;

    std::random_device                  random;
    std::mt19937                        mt(random());
    std::uniform_real_distribution<PT>  uniform(-100,100);

    auto func = [&](ET &val, Index i, Index j) mutable -> void
                {
                    val = ElementType(uniform(mt), uniform(mt));
                };

    hpc::matvec::apply(A, func);
}

//
//  C0 is the trusted result of C <- beta C + alpha*A*B
//  C1 is computed by a routine subject to testing
//
template <typename Alpha, typename MA, typename MB, typename Beta,
          typename MC0, typename MC1>
typename std::enable_if<hpc::matvec::IsGeMatrix<MA>::value
                     && hpc::matvec::IsGeMatrix<MB>::value
                     && hpc::matvec::IsGeMatrix<MC0>::value
                     && hpc::matvec::IsGeMatrix<MC1>::value,
         double>::type
estimateGemmResidual(const Alpha &alpha, const MA &A, const MB &B,
                     const Beta &beta, const MC0 &C0, const MC1 &C1)
{
    typedef typename MC0::ElementType   TC0;
    typedef typename MC0::Index         Index;

    Index m= C1.numRows;
    Index n= C1.numCols;
    Index k= A.numCols;

    double aNorm = hpc::matvec::asum(A) * std::abs(alpha);
    double bNorm = hpc::matvec::asum(B);
    double cNorm = hpc::matvec::asum(C1) * std::abs(beta);
    double diff = 0;
    hpc::matvec::apply(C0, [=,&diff](const TC0 &val, Index i, Index j) -> void
                           {
                               diff += std::abs(C1(i,j) - val);
                           });
    // Using eps for double gives upper bound in case elements have lower
    // precision.
    double eps = std::numeric_limits<double>::epsilon();
    double res = diff/(aNorm*bNorm*cNorm*eps*std::max(std::max(m,n),k));
    return res;
}

namespace refblas {

template <typename Index, typename Alpha, typename TA, typename TB,
          typename Beta, typename TC>
void
gemm(Index m, Index n, Index k,
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

} // namespace refblas

hpc::mt::ThreadPool global_tpool(4);

int
main()
{
    typedef double          Alpha;
    typedef double          TA;
    typedef double          TB;
    typedef double          Beta;
    typedef double          TC;
    typedef std::size_t     Index;

    hpc::matvec::StorageOrder order = (COLMAJOR)
                                    ? hpc::matvec::StorageOrder::ColMajor
                                    : hpc::matvec::StorageOrder::RowMajor;

    hpc::matvec::GeMatrix<TA, Index> A_(MAXDIM_M, MAXDIM_K, order);
    hpc::matvec::GeMatrix<TB, Index> B_(MAXDIM_K, MAXDIM_N, order);
    hpc::matvec::GeMatrix<TC, Index> C1_(MAXDIM_M, MAXDIM_N, order);
    hpc::matvec::GeMatrix<TC, Index> C2_(MAXDIM_M, MAXDIM_N, order);


    randomInit(A_);
    randomInit(B_);
    randomInit(C1_);

    const Alpha alpha(ALPHA);
    const Beta  beta(BETA);

    copy(C1_, C2_);

    // Header-Zeile fuer die Ausgabe
    fmt::printf("%5s %5s %5s ", "m", "n", "k");
    fmt::printf("%20s %9s", "refColMajor: t", "MFLOPS");
    fmt::printf("%20s %9s %9s", "blocked GEMM: t", "MFLOPS", "diff");
    fmt::printf("\n");

    hpc::aux::WallTime<double> wallTime;

    for (long m = MIN_M, n = MIN_N, k = MIN_K;
         m <=MAX_M && n <= MAX_N && k <= MAX_K;
         m += INC_M, n += INC_N, k += INC_K)
    {
        double t;

        auto A  = A_(0,0,m,k);
        auto B  = B_(0,0,k,n);
        auto C1 = C1_(0,0,m,n);
        auto C2 = C2_(0,0,m,n);

        fmt::printf("%5ld %5ld %5ld ", m, n, k);

        wallTime.tic();
        refblas::gemm(Index(m), Index(n), Index(k), alpha,
                      A.data, A.incRow, A.incCol,
                      B.data, B.incRow, B.incCol,
                      beta,
                      C1.data, C1.incRow, C1.incCol);
        t = wallTime.toc();
        fmt::printf("%20.4lf %9.2lf", t, 2.*m/1000*n/1000*k/t);

        wallTime.tic();
        hpc::matvec::mm(alpha, A, B, beta, C2);
        t = wallTime.toc();
        double res = estimateGemmResidual(alpha, A, B, beta, C1, C2);

        fmt::printf("%20.4lf %9.2lf %9.1e", t, 2.*m/1000*n/1000*k/t, res);
        fmt::printf("\n");
    }
}
