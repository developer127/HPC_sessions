#include <cassert>
#include <random>
#include <type_traits>
#include <hpc/matvec/densevector.hpp>
#include <hpc/matvec/iamax.hpp>
#include <hpc/matvec/r.hpp>
#include <hpc/matvec/scal.h>
#include <hpc/matvec/swap.hpp>
#include <hpc/matvec/print.h>
#include <hpc/ulmlapack/tf2.hpp>


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

template <typename VX>
typename std::enable_if<hpc::matvec::IsDenseVector<VX>::value,
                        void>::type
randomInit(VX &x)
{
    typedef typename VX::Index  Index;

    randomInit(x.length, Index(1), x.data, x.inc, Index(1));
}

template <typename MA>
typename std::enable_if<hpc::matvec::IsGeMatrix<MA>::value,
                        void>::type
randomInit(MA &A)
{
    randomInit(A.numRows, A.numCols, A.data, A.incRow, A.incCol);
}

//------------------------------------------------------------------------------

int
main()
{
    using namespace hpc::matvec;

    typedef double       T;
    typedef std::size_t  Size;
    typedef std::size_t  Index;

    GeMatrix<T, Size, Index>    A(8,17);
    auto x = A.row(2);
    auto y = x(3, 12);
    auto z = y(0, 4, 2);

    randomInit(A);

    print(A, "A");
    print(x, "x");
    print(y, "y");
    print(z, "z");

    printf("iamax(x) = %ld\n", iamax(x));
    printf("iamax(y) = %ld\n", iamax(y));
    printf("iamax(z) = %ld\n", iamax(z));

    auto a1 = A.row(0);
    scal(2.5, a1);

    auto x1 = x(0,5);
    auto x2 = x(5,5);

    swap(x1, x2);

    print(A, "A");
    print(x, "x");
    print(y, "y");
    print(z, "z");

    auto a10 = A.col(0)(1,7);
    auto a01 = A.row(0)(1,9);
    auto A11 = A(1,1,7,9);

    print(a01, "a01");
    print(a10, "a10");
    print(A11, "A11");

    r(1.5, a10, a01, A11);
    print(A, "A");

    auto p = A.col(0);

    //hpc::ulmlapack::tf2(A, p);
}
