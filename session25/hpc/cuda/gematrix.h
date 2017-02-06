#ifndef HPC_CUDA_GEMATRIX_H
#define HPC_CUDA_GEMATRIX_H 1

/* vim: set sw=4: */

#include <cassert>
#include <cstdlib>
#include <hpc/cuda/densevector.h>
#include <hpc/cuda/check.h>

namespace hpc { namespace cuda {

enum StorageOrder {
    ColMajor,
    RowMajor
};

template <typename T, typename I>
    struct GeMatrixView;

template <typename T, typename I>
    struct GeMatrixConstView;

template <typename T, typename I=std::size_t>
struct GeMatrix
{
    typedef T                              ElementType;
    typedef I                              Index;

    typedef GeMatrix<T,Index>              NoView;
    typedef GeMatrixConstView<T,Index>     ConstView;
    typedef GeMatrixView<T,Index>          View;

    typedef DenseVectorConstView<T,Index>  ConstVectorView;
    typedef DenseVectorView<T,Index>       VectorView;


    GeMatrix(Index numRows, Index numCols,
             StorageOrder order=ColMajor)
        : numRows(numRows), numCols(numCols),
          incRow(order==ColMajor ? 1: numCols),
          incCol(order==RowMajor ? 1: numRows)
    {
	CHECK_CUDA(cudaMalloc, (void**)&data, numRows*numCols * sizeof(T));
    }

    ~GeMatrix()
    {
	CHECK_CUDA(cudaFree, data);
    }

    ElementType*
    address(Index i, Index j)
    {
        assert(i<numRows && j<numCols);
        return &data[i*incRow + j*incCol];
    }

    const ElementType*
    address(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return &data[i*incRow + j*incCol];
    }

    ElementType &
    operator()(Index i, Index j)
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ConstView
    operator()(Index i, Index j, Index m, Index n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, address(i,j), incRow, incCol);
    }

    View
    operator()(Index i, Index j, Index m, Index n)
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return View(m, n, address(i,j), incRow, incCol);
    }

    ConstVectorView
    row(Index i) const
    {
        return ConstVectorView(numCols, address(i,0), incCol);
    }

    VectorView
    row(Index i)
    {
        return VectorView(numCols, address(i,0), incCol);
    }

    ConstVectorView
    col(Index j) const
    {
        return ConstVectorView(numRows, address(0,j), incRow);
    }

    VectorView
    col(Index j)
    {
        return VectorView(numRows, address(0,j), incRow);
    }

    const Index     numRows, numCols, incRow, incCol;
    ElementType*    data;

private:
    /* inhibit copy constructor */
    GeMatrix(const GeMatrix& other);
};

template <typename T, typename I=std::size_t>
struct GeMatrixView
{
    typedef T                               ElementType;
    typedef I                               Index;
    typedef GeMatrix<T,Index>               NoView;
    typedef GeMatrixConstView<T,Index>      ConstView;
    typedef GeMatrixView<T,Index>           View;

    typedef DenseVectorConstView<T,Index>   ConstVectorView;
    typedef DenseVectorView<T,Index>        VectorView;

    __device__ __host__
    GeMatrixView(Index numRows, Index numCols,
                 T *data,
                 Index incRow, Index incCol)
        : numRows(numRows), numCols(numCols),
          incRow(incRow), incCol(incCol), data(data)
    {
    }

    __device__ __host__
    GeMatrixView(const GeMatrixView &rhs)
        : numRows(rhs.numRows), numCols(rhs.numCols),
          incRow(rhs.incRow), incCol(rhs.incCol),
          data(rhs.data)
    {
    }

    __device__
    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    __device__
    ElementType &
    operator()(Index i, Index j)
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    __device__
    ConstView
    operator()(Index i, Index j, Index m, Index n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }


    __device__
    View
    operator()(Index i, Index j, Index m, Index n)
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return View(m, n, &(operator()(i,j)), incRow, incCol);
    }

    __device__
    ConstVectorView
    row(Index i) const
    {
        return ConstVectorView(numCols, &(operator()(i,0)), incCol);
    }

    __device__
    VectorView
    row(Index i)
    {
        return VectorView(numCols, &(operator()(i,0)), incCol);
    }

    __device__
    ConstVectorView
    col(Index j) const
    {
        return ConstVectorView(numRows, &(operator()(0,j)), incRow);
    }

    __device__
    VectorView
    col(Index j)
    {
        return VectorView(numRows, &(operator()(0,j)), incRow);
    }


    const Index     numRows, numCols, incRow, incCol;
    ElementType*    data;
};

template <typename T, typename I=std::size_t>
struct GeMatrixConstView
{
    typedef T                               ElementType;
    typedef I                               Index;
    typedef GeMatrix<T,Index>               NoView;
    typedef GeMatrixConstView<T,Index>      ConstView;
    typedef GeMatrixView<T,Index>           View;

    typedef DenseVectorConstView<T,Index>   ConstVectorView;
    typedef DenseVectorView<T,Index>        VectorView;

    __device__ __host__
    GeMatrixConstView(Index numRows, Index numCols,
                      const T *data,
                      Index incRow, Index incCol)
        : numRows(numRows), numCols(numCols),
          incRow(incRow), incCol(incCol), data(data)
    {
    }

    __device__ __host__
    GeMatrixConstView(const GeMatrixConstView &rhs)
        : numRows(rhs.numRows), numCols(rhs.numCols),
          incRow(rhs.incRow), incCol(rhs.incCol),
          data(rhs.data)
    {
    }

    __device__
    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    __device__
    ConstView
    operator()(Index i, Index j, Index m, Index n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }

    __device__
    ConstVectorView
    row(Index i) const
    {
        return ConstVectorView(numCols, &(operator()(i,0)), incCol);
    }

    __device__
    ConstVectorView
    col(Index j) const
    {
        return ConstVectorView(numRows, &(operator()(0,j)), incRow);
    }

    const Index         numRows, numCols, incRow, incCol;
    const ElementType*  data;
};



} } // namespace cuda, hpc

#endif // HPC_CUDA_GEMATRIX_H
