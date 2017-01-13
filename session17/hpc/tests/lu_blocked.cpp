#include <cassert>
#include <random>
#include <fmt/printf.hpp>
#include <hpc/aux/walltime.h>
#include <hpc/matvec/copy.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/mm.h>
#include <hpc/matvec/print.h>
#include <hpc/ulmlapack/getrf.hpp>

using namespace hpc::matvec;

template <typename T, typename I>
double
checkLU(const GeMatrixView<T, I> &A,
        GeMatrixView<T, I> &A_,
        const GeMatrixView<I, I> &p)
{
    assert(A.numRows==A_.numRows);
    assert(A.numCols==A_.numCols);
    assert(A.numCols>=A.numRows);
    assert(A.numRows==p.numRows);

    typedef typename GeMatrix<T, I>::Index  Index;

    Index m = A.numRows;
    Index n = A.numRows;

    // Copy L form A_ and set A_ to U
    GeMatrix<T, I> L(m, m);
    for (Index j=0; j<m; ++j) {
        for (Index i=0; i<m; ++i) {
            if (i==j) {
                L(i,i) = 1;
            } else if (i>j) {
                L(i,j)  = A_(i,j);
                A_(i,j) = 0;
            } else if (j>i) {
                L(i,j)  = 0;
            }
        }
    }

    // Compute LU = L*U
    GeMatrix<T, I> LU(m, n);
    for (Index j=0; j<n; ++j) {
        for (Index i=0; i<m; ++i) {
            LU(i,j) = 0;
        }
    }
    hpc::matvec::mm(T(1), L, A_, T(0), LU);

    // Apply P
    for (Index i=m; i>0; --i) {
        if (i-1!=p(i-1,0)) {
            hpc::ulmblas::swap(n,
                               &LU(i-1,0), LU.incCol,
                               &LU(p(i-1,0),0), LU.incCol);
        }
    }

    double diff = 0;
    for (Index j=0; j<n; ++j) {
        for (Index i=0; i<m; ++i) {
            diff += std::abs(A(i,j)-LU(i,j));
        }
    }
    return diff;
}

//
//  Random initializer for general matrices: real and complex valued
//
template <typename Index, typename T>
void
randomInit(Index m, Index n, T *A, Index incRowA, Index incColA)
{
    std::random_device                  random;
    std::default_random_engine          mt(random());
    std::uniform_real_distribution<T>   uniform(-100,100);

    for (Index i=0; i<m; ++i) {
        for (Index j=0; j<n; ++j) {
            A[i*incRowA+j*incColA] = uniform(mt);
        }
    }
}

int
main()
{
    const std::ptrdiff_t   m_min = 10;
    const std::ptrdiff_t   n_min = 10;

    const std::ptrdiff_t   m_max = 2000;
    const std::ptrdiff_t   n_max = 2000;

    const std::ptrdiff_t   m_inc = 10;
    const std::ptrdiff_t   n_inc = 10;

    GeMatrix<double, std::ptrdiff_t>      A(m_max, n_max);
    GeMatrix<double, std::ptrdiff_t>      A_org(m_max, n_max);
    GeMatrix<std::ptrdiff_t, std::ptrdiff_t>        p(m_max, 1);
    hpc::aux::WallTime<double>  timer;

    for (std::ptrdiff_t m=m_min, n=n_min; m<=m_max && n<=n_max; m+=m_inc, n+=n_inc) {

        randomInit(m, n, A.data, A.incRow, A.incCol);
        hpc::ulmblas::gecopy(m, n,
                             A.data, A.incRow, A.incCol,
                             A_org.data, A_org.incRow, A_org.incCol);
        //print(A_org, "A_org");

        timer.tic();
        std::ptrdiff_t info = hpc::ulmlapack::getrf(m, n,
                                                    A.data, A.incRow, A.incCol,
                                                    p.data, p.incRow);
        double t = timer.toc();
        //print(A, "A");
        //print(p, "p");

        auto A_ = A_org(0, 0, m, n);
        auto LU = A(0, 0, m, n);
        auto p_ = p(0, 0, m, 1);

        double diff = checkLU(A_, LU, p_);
        //print(A, "LU");
        fmt::printf("%4ld %4ld %4ld %16.3e %16.5lf\n", m, n, info, diff, t);
    }
}
