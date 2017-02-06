#ifndef HPC_CUDA_DENSEVECTOR_H
#define HPC_CUDA_DENSEVECTOR_H 1

/* vim: set sw=4: */

#ifndef __CUDACC__
#	error This source must be compiled using nvcc
#endif

#include <cassert>
#include <cstdlib>
#include <hpc/cuda/check.h>

namespace hpc { namespace cuda {

template <typename T, typename I>
    struct DenseVectorConstView;

template <typename T, typename I>
    struct DenseVectorView;

template <typename T, typename I=std::size_t>
struct DenseVector
{
    typedef T                               ElementType;
    typedef I                               Index;

    typedef DenseVector<T,Index>            NoView;
    typedef DenseVectorConstView<T,Index>   ConstView;
    typedef DenseVectorView<T,Index>        View;

    DenseVector(Index length)
        : length(length), inc(1)
    {
	CHECK_CUDA(cudaMalloc, (void**)&data, length * sizeof(T));
    }

    ~DenseVector()
    {
	CHECK_CUDA(cudaFree, data);
    }

    ElementType*
    address(Index i)
    {
	return &data[i*inc];
    }

    const ElementType*
    address(Index i) const
    {
	return &data[i*inc];
    }

    ConstView
    operator()(Index i, Index num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, address(i), stride);
    }

    View
    operator()(Index i, Index num, Index stride=1)
    {
        assert(i+stride*num<=length);
        return View(num, address(i), stride);
    }

    const Index     length, inc;
    ElementType*    data;

private:
    /* inhibit copy constructor */
    DenseVector(const DenseVector& other);
};

template <typename T, typename I=std::size_t>
struct DenseVectorView
{
    typedef T                           ElementType;
    typedef I                           Index;

    typedef DenseVector<T,Index>            NoView;
    typedef DenseVectorConstView<T,Index>   ConstView;
    typedef DenseVectorView<T,Index>        View;

    __device__ __host__
    DenseVectorView(Index length, ElementType *data, Index inc)
        : length(length), inc(inc), data(data)
    {
    }

    __device__
    const ElementType &
    operator()(Index i) const
    {
        assert(i<length);
        return data[i*inc];
    }

    __device__
    ElementType &
    operator()(Index i)
    {
        assert(i<length);
        return data[i*inc];
    }

    __device__
    ConstView
    operator()(Index i, Index num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, &(operator()(i)), inc*stride);
    }

    __device__
    View
    operator()(Index i, Index num, Index stride=1)
    {
        assert(i+stride*num<=length);
        return View(num, &(operator()(i)), inc*stride);
    }

    const Index     length, inc;
    ElementType*    data;
};

template <typename T, typename I=std::size_t>
struct DenseVectorConstView
{
    typedef T                           ElementType;
    typedef I                           Index;

    typedef DenseVector<T,Index>            NoView;
    typedef DenseVectorConstView<T,Index>   ConstView;
    typedef DenseVectorView<T,Index>        View;

    __device__ __host__
    DenseVectorConstView(Index length, const ElementType *data, Index inc)
        : length(length), inc(inc), data(data)
    {
    }

    __device__
    const ElementType &
    operator()(Index i) const
    {
        assert(i<length);
        return data[i*inc];
    }

    __device__
    ConstView
    operator()(Index i, Index num, Index stride=1) const
    {
        assert(i+stride*num<=length);
        return ConstView(num, &(operator()(i)), inc*stride);
    }

    const Index         length, inc;
    const ElementType*  data;
};


} } // namespace cuda, hpc

#endif // HPC_CUDA_DENSEVECTOR_H
