#include <cassert>
#include <random>
#include <type_traits>
#include <hpc/matvec/densevector.hpp>
#include <hpc/matvec/isdensevector.hpp>
#include <hpc/matvec/print.h>


//
//  Random initializer for general matrices: real and complex valued
//
template <typename Size, typename Index, typename T>
void
randomInit(Size m, Size n, T *A, Index incRowA, Index incColA)
{
    std::random_device                  random;
    std::default_random_engine          mt(random());
    std::uniform_real_distribution<T>   uniform(-100,100);

    for (Size i=0; i<m; ++i) {
        for (Size j=0; j<n; ++j) {
            A[i*incRowA+j*incColA] = uniform(mt);
        }
    }
}

template <typename VX>
typename std::enable_if<hpc::matvec::IsDenseVector<VX>::value,
                        void>::type
randomInit(VX &x)
{
    typedef typename VX::Size   Size;
    typedef typename VX::Index  Index;

    randomInit(x.length, Size(1), x.data, x.inc, Index(1));
}

//------------------------------------------------------------------------------

template <typename VX>
typename std::enable_if<hpc::matvec::IsDenseVector<VX>::value,
                        void>::type
foo(const VX &x)
{
    std::printf("Entering foo\n");
    auto y = x(1, 12);
    auto z = y(0, 3, 2);

    print(x, "x");
    print(y, "y");
    print(z, "z");
    std::printf("Leaving foo\n");
}

int
main()
{
    using namespace hpc::matvec;

    typedef double       T;
    typedef std::size_t  Index;

    DenseVector<T, Index> x(20);

    auto y = x(1, 12);
    auto z = y(0, 3, 2);

    randomInit(x);

    print(x, "x");
    print(y, "y");
    print(z, "z");

    foo(x);
}
