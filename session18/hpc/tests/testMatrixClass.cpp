#include <cassert>
#include <random>
#include <type_traits>
#include <hpc/matvec/densevector.hpp>
#include <hpc/matvec/print.h>


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

template <typename MA>
typename std::enable_if<hpc::matvec::IsGeMatrix<MA>::value,
                        void>::type
foo(const MA &A)
{
    std::printf("Entering foo\n");

    auto m = A.numRows;
    auto n = A.numCols;

    auto B = A(0, 0, m, n);

    auto x = B.row(0);
    auto y = B.col(0);

    print(x, "x");
    print(y, "y");
    std::printf("Leaving foo\n");
}

int
main()
{
    using namespace hpc::matvec;

    typedef double       T;
    typedef std::size_t  Index;

    GeMatrix<T, Index> A(8,10);

    auto x = A.row(3);
    auto y = A.col(4);

    randomInit(A);

    print(A, "A");
    print(x, "x");
    print(y, "y");

    auto B = A(1, 2, 6, 5);

    auto v = B.row(3);
    auto w = B.col(4);

    print(B, "B");
    print(v, "v");
    print(w, "w");

    foo(B);
}
