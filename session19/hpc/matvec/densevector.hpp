#ifndef HPC_MATVEC_DENSEVECTOR_HPP
#define HPC_MATVEC_DENSEVECTOR_HPP 1

#include <cassert>
#include <cstdlib>
#include <hpc/matvec/copy.h>

namespace hpc { namespace matvec {

template <typename T, typename S, typename I>
    struct DenseVectorConstView;

template <typename T, typename S, typename I>
    struct DenseVectorView;

template <typename T, typename S=std::size_t, typename I=std::ptrdiff_t>
struct DenseVector
{
    typedef T                                   ElementType;
    typedef S                                   Size;
    typedef I                                   Index;

    typedef DenseVector<T,Size,Index>           NoView;
    typedef DenseVectorConstView<T,Size,Index>  ConstView;
    typedef DenseVectorView<T,Size,Index>       View;

    DenseVector(Size length): length(length), inc(1), data(new T[length])
    {
    }

    ~DenseVector()
    {
        delete[] data;
    }

    template <typename VX>
    typename std::enable_if<IsDenseVector<VX>::value,
        NoView>::type &
    operator=(VX &rhs)
    {
        if(this != &rhs){
            copy(rhs, *this);
        }
        return *this;
    }


    const ElementType &
    operator()(Size i) const
    {
        assert(i<length);
        return data[i*inc];
    }

    ElementType &
    operator()(Size i)
    {
        assert(i<length);
        return data[i*inc];
    }

    ConstView
    operator()(Size i, Size num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, &operator()(i), inc*stride);
    }

    View
    operator()(Size i, Size num, Index stride=1)
    {
        assert(i+stride*num<=length);
        return View(num, &operator()(i), inc*stride);
    }

    const Size      length;
    const Index     inc;
    ElementType*    data;

};

template <typename T, typename S=std::size_t, typename I=std::ptrdiff_t>
struct DenseVectorView
{
    typedef T                                   ElementType;
    typedef S                                   Size;
    typedef I                                   Index;

    typedef DenseVector<T,Size,Index>           NoView;
    typedef DenseVectorConstView<T,Size,Index>  ConstView;
    typedef DenseVectorView<T,Size,Index>       View;

    DenseVectorView(Size length, ElementType *data, Index inc)
        : length(length), inc(inc), data(data)
    {
    }

    DenseVectorView(const DenseVectorView &rhs)
        : length(rhs.length), inc(rhs.inc), data(rhs.data)
    {
    }

    const ElementType &
    operator()(Size i) const
    {
        assert(i<length);
        return data[i*inc];
    }

    ElementType &
    operator()(Size i)
    {
        assert(i<length);
        return data[i*inc];
    }

    ConstView
    operator()(Size i, Size num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, &operator()(i), inc*stride);
    }

    View
    operator()(Size i, Size num, Index stride=1)
    {
        assert(i+stride*num<=length);
        return View(num, &operator()(i), inc*stride);
    }

    const Size      length;
    const Index     inc;
    ElementType*    data;
};

template <typename T, typename S=std::size_t, typename I=std::ptrdiff_t>
struct DenseVectorConstView
{
    typedef T                                   ElementType;
    typedef S                                   Size;
    typedef I                                   Index;

    typedef DenseVector<T,Size,Index>           NoView;
    typedef DenseVectorConstView<T,Size,Index>  ConstView;
    typedef DenseVectorView<T,Size,Index>       View;

    DenseVectorConstView(Size length, const ElementType *data, Index inc)
        : length(length), inc(inc), data(data)
    {
    }

    DenseVectorConstView(const DenseVectorConstView &rhs)
        : length(rhs.length), inc(rhs.inc), data(rhs.data)
    {
    }

    const ElementType &
    operator()(Size i) const
    {
        assert(i<length);
        return data[i*inc];
    }

    ConstView
    operator()(Size i, Size num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, &operator()(i), inc*stride);
    }

    const Size      length;
    const Index     inc;
    const ElementType*  data;
};


} } // namespace matvec, hpc

#endif // HPC_MATVEC_DENSEVECTOR_HPP
