#include <cassert>
#include <random>
#include <printf.hpp>
#include <hpc/matvec/copy.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/mm.h>
#include <hpc/matvec/print.h>
#include <hpc/ulmblas/trlsm.hpp>
#include <hpc/ulmblas/trusm.hpp>
#include <hpc/ulmlapack/getf2.hpp>
#include <hpc/ulmlapack/swap.hpp>


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
    using namespace hpc::matvec;

    typedef double       T;
    typedef std::size_t  Index;

    const Index   m = 8;
    const Index   n = 8;

    GeMatrix<T>      A(m, m);
    GeMatrix<T>      X(m, n);
    GeMatrix<T>      B(m, n);
    GeMatrix<Index>  p(m, 1);

    randomInit(m, m, A.data, A.incRow, A.incCol);
    randomInit(m, n, X.data, X.incRow, X.incCol);

    //
    //  A*X = B
    //
    mm(T(1), A, X, T(0), B);

    //
    //  Factorize A=P*L*U
    //
    std::ptrdiff_t info = hpc::ulmlapack::getf2(A.numRows, A.numCols,
                                             A.data, A.incRow, A.incCol,
                                             p.data, p.incRow);
    fmt::printf("getf2 returned: info = %ld\n", info);

    //
    //  Solve P*L*U*X = B
    //

    print(p, "Pivot");
    hpc::ulmlapack::swap(B.numRows, B.numCols, B.data, B.incRow, B.incCol,
                         Index(0), B.numRows,
                         p.data, p.incRow);
    hpc::ulmblas::trlsm(m, n, T(1), true,
                        A.data, A.incRow, A.incCol,
                        B.data, B.incRow, B.incCol);
    hpc::ulmblas::trusm(m, n, T(1), false,
                        A.data, A.incRow, A.incCol,
                        B.data, B.incRow, B.incCol);


    //
    //  Compute residual
    //
    double res = 0;
    for (Index j=0; j<n; ++j) {
        for (Index i=0; i<m; ++i) {
            res += std::abs(B(i,j)-X(i,j));
        }
    }

    fmt::printf("res=%e\n", res);
}
