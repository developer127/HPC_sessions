#ifndef HPC_AUX_PRIMITIVE_TYPE_H
#define HPC_AUX_PRIMITIVE_TYPE_H 1

#include <complex>

namespace hpc { namespace aux {

template <typename T>
struct PrimitiveType
{
    typedef T       Type;
};

template <typename T>
struct PrimitiveType<std::complex<T> >
{
    typedef T       Type;
};

} } // namespace aux, hpc

#endif // HPC_AUX_PRIMITIVE_TYPE_H
