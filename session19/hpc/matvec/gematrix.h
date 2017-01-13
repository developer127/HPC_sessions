#ifndef HPC_MATVEC_GEMATRIX_H
#define HPC_MATVEC_GEMATRIX_H 1

#include <cassert>
#include <hpc/matvec/densevector.hpp>

namespace hpc { namespace matvec {

enum class StorageOrder {
    ColMajor,
    RowMajor
};

template <typename T, typename S, typename I>
    struct GeMatrixView;

template <typename T, typename S, typename I>
    struct GeMatrixConstView;

template <typename T, typename S, typename I>
    struct DenseVectorConstView;

template <typename T, typename S, typename I>
    struct DenseVectorView;

template <typename T, typename S=std::size_t, typename I=std::size_t>
struct GeMatrix
{
    typedef T                                   ElementType;
    typedef S                                   Size;
    typedef I                                   Index;
    typedef GeMatrix<T,Size,Index>              NoView;
    typedef GeMatrixConstView<T,Size,Index>     ConstView;
    typedef GeMatrixView<T,Size,Index>          View;

    typedef DenseVectorConstView<T,Size,Index>  ConstVectorView;
    typedef DenseVectorView<T,Size,Index>       VectorView;

    GeMatrix(Size numRows, Size numCols,
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
    operator()(Size i, Size j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ElementType &
    operator()(Size i, Size j)
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ConstView
    operator()(Size i, Size j, Size m, Size n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }

    View
    operator()(Size i, Size j, Size m, Size n)
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return View(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Size i) const
    {
        assert(i<numRows);
        return ConstVectorView(numCols, &data[i*incRow], incCol);
    }

    VectorView
    row(Size i)
    {
        assert(i<numRows);
        return VectorView(numCols, &data[i*incRow], incCol);
    }

    ConstVectorView
    col(Size j) const
    {
        assert(j<numRows);
        return ConstVectorView(numRows, &data[j*incCol], incRow);
    }

    VectorView
    col(Size j)
    {
        assert(j<numRows);
        return VectorView(numRows, &data[j*incCol], incRow);
    }

    const Size      numRows, numCols;
    const Index     incRow, incCol;
    ElementType*    data;
};

template <typename T, typename S=std::size_t, typename I=std::size_t>
struct GeMatrixView
{
    typedef T                                   ElementType;
    typedef S                                   Size;
    typedef I                                   Index;
    typedef GeMatrix<T,Size,Index>              NoView;
    typedef GeMatrixConstView<T,Size,Index>     ConstView;
    typedef GeMatrixView<T,Size,Index>          View;

    typedef DenseVectorConstView<T,Size,Index>  ConstVectorView;
    typedef DenseVectorView<T,Size,Index>       VectorView;

    GeMatrixView(Size numRows, Size numCols,
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
    operator()(Size i, Size j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ElementType &
    operator()(Size i, Size j)
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ConstView
    operator()(Size i, Size j, Size m, Size n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }


    View
    operator()(Size i, Size j, Size m, Size n)
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return View(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Size i) const
    {
        assert(i<numRows);
        return ConstVectorView(numCols, &data[i*incRow], incCol);
    }

    VectorView
    row(Size i)
    {
        assert(i<numRows);
        return VectorView(numCols, &data[i*incRow], incCol);
    }

    ConstVectorView
    col(Size j) const
    {
        assert(j<numRows);
        return ConstVectorView(numRows, &data[j*incCol], incRow);
    }

    VectorView
    col(Size j)
    {
        assert(j<numRows);
        return VectorView(numRows, &data[j*incCol], incRow);
    }

    const Size      numRows, numCols;
    const Index     incRow, incCol;
    ElementType*    data;
};

template <typename T, typename S=std::size_t, typename I=std::size_t>
struct GeMatrixConstView
{
    typedef T                                   ElementType;
    typedef S                                   Size;
    typedef I                                   Index;
    typedef GeMatrix<T,Size,Index>              NoView;
    typedef GeMatrixConstView<T,Size,Index>     ConstView;
    typedef GeMatrixView<T,Size,Index>          View;

    typedef DenseVectorConstView<T,Size,Index>  ConstVectorView;

    GeMatrixConstView(Size numRows, Size numCols,
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
    operator()(Size i, Size j) const
    {
        assert(i<numRows && j<numCols);
        return data[i*incRow + j*incCol];
    }

    ConstView
    operator()(Size i, Size j, Size m, Size n) const
    {
        assert(i+m<=numRows);
        assert(j+n<=numCols);
        return ConstView(m, n, &(operator()(i,j)), incRow, incCol);
    }

    ConstVectorView
    row(Size i) const
    {
        assert(i<numRows);
        return ConstVectorView(numCols, &data[i*incRow], incCol);
    }

    ConstVectorView
    col(Size j) const
    {
        assert(j<numRows);
        return ConstVectorView(numRows, &data[j*incCol], incRow);
    }

    const Size          numRows, numCols;
    const Index         incRow, incCol;
    const ElementType*  data;
};



} } // namespace matvec, hpc

#endif // HPC_MATVEC_GEMATRIX_H
