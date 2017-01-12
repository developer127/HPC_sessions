#ifndef HPC_MATVEC_GEMATRIX_H
#define HPC_MATVEC_GEMATRIX_H 1

#include <cassert>
#include <hpc/matvec/densevector.hpp>

namespace hpc { namespace matvec {

enum class StorageOrder {
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
    typedef T                           ElementType;
    typedef I                           Index;
    typedef GeMatrix<T,Index>           NoView;
    typedef GeMatrixConstView<T,Index>  ConstView;
    typedef GeMatrixView<T,Index>       View;

    typedef DenseVectorView<T,std::size_t,Index>       DenseVectorView;

    GeMatrix(Index numRows, Index numCols,
             StorageOrder order=StorageOrder::ColMajor)
        : numRows(numRows), numCols(numCols),
          incRow(order==StorageOrder::ColMajor ? 1: numCols),
          incCol(order==StorageOrder::RowMajor ? 1: numRows),
          data(new T[numRows*numCols])
    {
    }

    ~GeMatrix()
    {
        delete[] data;
    }

    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
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
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }

    View
    operator()(Index i, Index j, Index m, Index n)
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return View(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Index i) const
    {
        assert(i<numRows);
        return DenseVectorView(numCols, &data[i*incRow], incCol);
    }

    DenseVectorView
    row(Index i)
    {
        assert(i<numRows);
        return DenseVectorView(numCols, &data[i*incRow], incCol);
    }

    ConstVectorView
    col(Index j) const
    {
        assert(i<numRows);
        return DenseVectorView(numRows, &data[j*incCol], incRow);
    }

    VectorView
    col(Index j)
    {
        assert(i<numRows);
        return DenseVectorView(numRows, &data[j*incCol], incRow);
    }

    const Index     numRows, numCols, incRow, incCol;
    ElementType*    data;
};

template <typename T, typename I=std::size_t>
struct GeMatrixView
{
    typedef T                           ElementType;
    typedef I                           Index;
    typedef GeMatrix<T,Index>           NoView;
    typedef GeMatrixConstView<T,Index>  ConstView;
    typedef GeMatrixView<T,Index>       View;

    GeMatrixView(Index numRows, Index numCols,
                 T *data,
                 Index incRow, Index incCol)
        : numRows(numRows), numCols(numCols),
          incRow(incRow), incCol(incCol), data(data)
    {
    }

    GeMatrixView(const GeMatrixView &rhs)
        : numRows(rhs.numRows), numCols(rhs.numCols),
          incRow(rhs.incRow), incCol(rhs.incCol),
          data(rhs.data)
    {
    }

    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
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
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }


    View
    operator()(Index i, Index j, Index m, Index n)
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return View(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Index i) const
    {
        assert(i<numRows);
        return DenseVectorView(numCols, &data[i*incRow], incCol);
    }

    DenseVectorView
    row(Index i)
    {
        assert(i<numRows);
        return DenseVectorView(numCols, &data[i*incRow], incCol);
    }

    ConstVectorView
    col(Index j) const
    {
        assert(i<numRows);
        return DenseVectorView(numRows, &data[j*incCol], incRow);
    }

    VectorView
    col(Index j)
    {
        assert(i<numRows);
        return DenseVectorView(numRows, &data[j*incCol], incRow);
    }

    const Index     numRows, numCols, incRow, incCol;
    ElementType*    data;
};

template <typename T, typename I=std::size_t>
struct GeMatrixConstView
{
    typedef T                           ElementType;
    typedef I                           Index;
    typedef GeMatrix<T,Index>           NoView;
    typedef GeMatrixConstView<T,Index>  ConstView;
    typedef GeMatrixView<T,Index>       View;

    GeMatrixConstView(Index numRows, Index numCols,
                      const T *data,
                      Index incRow, Index incCol)
        : numRows(numRows), numCols(numCols),
          incRow(incRow), incCol(incCol), data(data)
    {
    }

    GeMatrixConstView(const GeMatrixConstView &rhs)
        : numRows(rhs.numRows), numCols(rhs.numCols),
          incRow(rhs.incRow), incCol(rhs.incCol),
          data(rhs.data)
    {
    }

    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ConstView
    operator()(Index i, Index j, Index m, Index n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Index i) const
    {
        assert(i<numRows);
        return DenseVectorView(numCols, &data[i*incRow], incCol);
    }

    ConstVectorView
    col(Index j) const
    {
        assert(i<numRows);
        return DenseVectorView(numRows, &data[j*incCol], incRow);
    }

    const Index         numRows, numCols, incRow, incCol;
    const ElementType*  data;
};



} } // namespace matvec, hpc

#endif // HPC_MATVEC_GEMATRIX_H
