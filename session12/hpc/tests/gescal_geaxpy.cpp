#include <cstdio>
#include <complex>
#include <hpc/ulmblas/geaxpy.hpp>
#include <hpc/ulmblas/gecopy.hpp>
#include <hpc/ulmblas/gescal.hpp>

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
    printf(" %6.1f", x);
}

void
print_value(double x)
{
    printf(" %6.1lf", x);
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
    // Element type for A, B and alpha
    typedef float   TA;
    typedef double  TB;
    typedef double  Alpha;

    // Dimensions of A and B 
    std::size_t m       = 5;
    std::size_t n       = 6;

    // Storage order of A
    std::ptrdiff_t incRowA = 1;
    std::ptrdiff_t incColA = m;

    // Storage order of B
    std::ptrdiff_t incRowB = 1;
    std::ptrdiff_t incColB = m;

    // Allocate A and B
    TA *A = new TA[m*n];
    TB *B = new TB[m*n];

    // Value for alpha
    Alpha alpha = 1.5;
    printf("alpha =\n");
    geprint(1, 1, &alpha, 1, 1);

    // Initialize A
    auto fooInit = [](std::size_t m, std::size_t n,
                      std::size_t i, std::size_t j) -> TA
                   {
                       return TA(m)*TA(j)+TA(i)+TA(1);
                   };
    geinit(m, n, A, incRowA, incColA, fooInit);

    // Print A
    printf("A = \n");
    geprint(m, n, A, incRowA, incColA);

    // Copy A to B
    hpc::ulmblas::gecopy(m, n, A, incRowA, incColA, B, incRowB, incColB);
    printf("B = \n");
    geprint(m, n, B, incRowB, incColB);

    //  B <- alpha*B
    printf("B <- alpha*B\n");
    printf("B = \n");
    hpc::ulmblas::gescal(m, n, alpha, B, incRowB, incColB);
    geprint(m, n, B, incRowB, incColB);

    //  B <- alpha*A + B
    printf("B <- alpha*A +B\n");
    printf("B = \n");
    hpc::ulmblas::geaxpy(m, n, alpha, A, incRowA, incColA, B, incRowB, incColB);
    geprint(m, n, B, incRowB, incColB);
}
