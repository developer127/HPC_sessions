#ifndef HPC_MATVEC_PRINT_H
#define HPC_MATVEC_PRINT_H 1

#include <complex>
#include <type_traits>
#include <hpc/matvec/isgematrix.h>
#include <printf.hpp> // from https://github.com/afborchert/fmt

namespace hpc { namespace matvec {

template <typename T>
struct Format
{
    static const char *
    value()
    {
	// as fmt::printf is type-safe this can be used as catch-all
	// for all types
        return " %10.1f";
    }
};

template <>
struct Format<std::complex<double> >
{
    static const char *
    value()
    {
        return " (%10.1lf, %10.1lf)";
    }
};

//------------------------------------------------------------------------------

template <typename T>
void
print_value(T value) {
   fmt::printf(Format<T>::value(), value);
}

template <typename T>
void
print_value(const std::complex<T> &value) {
   fmt::printf(Format<std::complex<T> >::value(), value.real(), value.imag());
}

template <typename MA>
typename std::enable_if<IsGeMatrix<MA>::value,
         void>::type
print(const MA &A, const char *name = "")
{
    typedef typename MA::Index    Index;

    if (*name) {
        fmt::printf("%s = \n", name);
    }
    for (Index i=0; i<A.numRows; ++i) {
        for (Index j=0; j<A.numCols; ++j) {
            print_value(A(i,j));
        }
        fmt::printf("\n");
    }
    fmt::printf("\n");
}

} } // namespace matvec, hpc

#endif // HPC_MATVEC_PRINT_H
