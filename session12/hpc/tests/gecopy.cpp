#include <cstdio>
#include <complex>
#include <hpc/ulmblas/gecopy.hpp>

template <typename T, typename Func, typename Size, typename Index>
void
geinit(Size m, Size n, T *A, Index incRowA, Index incColA, Func func)
{
    for (Size i=0; i<m; ++i) {
        for (Size j=0; j<n; ++j) {
            A[i*incRowA+j*incColA] = func(m, n, i, j);
        }
    }
}

void
print_value(float x)
{
    printf("%6.1f", x);
}

void
print_value(double x)
{
    printf("%6.1lf", x);
}

template <typename T>
void
print_value(std::complex<T> z)
{
    printf("(");
    print_value(z.real());
    printf(",");
    print_value(z.imag());
    printf(")");
}

template <typename T, typename Size, typename Index>
void
geprint(Size m, Size n, const T *A, Index incRowA, Index incColA)
{
    for (Size i=0; i<m; ++i) {
        for (Size j=0; j<n; ++j) {
            print_value(A[i*incRowA+j*incColA]);
        }
        printf("\n");
    }
}

int
main()
{
    // Element type for A and B
    typedef float  TB;
    typedef double  TA;

    // Dimensions of A and B 
    std::size_t m       = 7;
    std::size_t n       = 8;

    // Storage order of A
    std::ptrdiff_t incRowA = 1;
    std::ptrdiff_t incColA = m;

    // Storage order of B
    std::ptrdiff_t incRowB = n;
    std::ptrdiff_t incColB = 1;

    // Allocate A and B
    TA *A = new TA[m*n];
    TB *B = new TB[m*n];

    // Initialize A
    auto fooInit = [](std::size_t m, std::size_t n,
                      std::size_t i, std::size_t j) -> TA
                   {
                       return TA(m)*TA(j)+TA(i)+TA(1);
                   };
    geinit(m, n, A, incRowA, incColA, fooInit);

    // Copy A to B
    hpc::ulmblas::gecopy(m, n, A, incRowA, incColA, B, incRowB, incColB);

    // Print B
    geprint(m, n, B, incRowB, incColB);
}
