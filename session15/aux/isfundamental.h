#ifndef HPC_AUX_ISFUNDAMENTAL_H
#define HPC_AUX_ISFUNDAMENTAL_H 1

#include <complex>
#include <type_traits>

namespace hpc { namespace aux {

template <typename T>
struct IsFundamental
{
    static const bool value = std::is_same<T, int>::value
                           || std::is_same<T, long>::value
                           || std::is_same<T, float>::value
                           || std::is_same<T, double>::value
                           || std::is_same<T, std::complex<float> >::value
                           || std::is_same<T, std::complex<double> >::value;
};

} } // namespace aux, hpc

#endif // HPC_AUX_ISFUNDAMENTAL_H
