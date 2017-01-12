#ifndef HPC_MATVEC_ISDENSEVECTOR_HPP
#define HPC_MATVEC_ISDENSEVECTOR_HPP 1

#include <cassert>
#include <type_traits>
#include <hpc/aux/iscomplex.h>
#include <hpc/matvec/densevector.hpp>

namespace hpc { namespace matvec {

template <typename Any>
struct IsDenseVector_
{
    static constexpr bool value = false;
};

template <typename T, typename S, typename I>
struct IsDenseVector_<DenseVector<T,S,I> >
{
    static constexpr bool value = true;
};

template <typename T, typename S, typename I>
struct IsDenseVector_<DenseVectorView<T,S,I> >
{
    static constexpr bool value = true;
};

template <typename T, typename S, typename I>
struct IsDenseVector_<DenseVectorConstView<T,S,I> >
{
    static constexpr bool value = true;
};

template <typename Any_>
struct IsDenseVector
{
    typedef typename std::remove_reference<Any_>::type  Any;
    static constexpr bool value = IsDenseVector_<Any>::value;
};

template <typename Vector>
struct IsRealDenseVector
{
    typedef typename Vector::ElementType    T;

    static constexpr bool value = IsDenseVector<Vector>::value
                               && !aux::IsComplex<T>::value;
};

template <typename Vector>
struct IsComplexDenseVector
{
    typedef typename Vector::ElementType    T;

    static constexpr bool value = IsDenseVector<Vector>::value
                               && aux::IsComplex<T>::value;
};

} } // namespace matvec, hpc

#endif // HPC_MATVEC_ISDENSEVECTOR_HPP
