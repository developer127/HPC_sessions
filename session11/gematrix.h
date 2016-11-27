#ifndef HPC_GEMATRIX_H
#define HPC_GEMATRIX_H 1

#include <cstddef>
#include <cassert>
#include "bench.h"
#include <functional>

namespace matvec {

enum class StorageOrder {
    ColMajor,
    RowMajor
};

template <typename T, typename I>
    struct GeMatrixView;

template <typename T, typename I=std::size_t>
struct GeMatrix
{
    typedef T                       ElementType;
    typedef I                       Index;
    typedef GeMatrix<T,Index>       NoView;
    typedef GeMatrixView<T,Index>   View;

    GeMatrix(Index m, Index n, StorageOrder order=StorageOrder::ColMajor)
        : m(m), n(n),
          incRow(order==StorageOrder::ColMajor ? 1: n),
          incCol(order==StorageOrder::RowMajor ? 1: m),
          data(new T[m*n])
    {
    }

    ~GeMatrix()
    {
        delete[] data;
    }

    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<m && j<n);
        return data[i*incRow + j*incCol];
    }

    ElementType &
    operator()(Index i, Index j)
    {
        assert(i<m && j<n);
        return data[i*incRow + j*incCol];
    }

    View
    operator()(Index i, Index j, Index numRows, Index numCols)
    {
        assert(i+numRows<=m);
        assert(j+numCols<=n);
        return View(numRows, numCols, &(operator()(i,j)), incRow, incCol);
    }

    const Index     m, n, incRow, incCol;
    ElementType*    data;
};

template <typename T, typename I=std::size_t>
struct GeMatrixView
{
    typedef T                       ElementType;
    typedef I                       Index;
    typedef GeMatrix<T,Index>       NoView;
    typedef GeMatrixView<T,Index>   View;

    GeMatrixView(Index m, Index n, T *data, Index incRow, Index incCol)
        : m(m), n(n), incRow(incRow), incCol(incCol), data(data)
    {
    }

    GeMatrixView(const GeMatrixView &rhs)
        : m(rhs.m), n(rhs.n), incRow(rhs.incRow), incCol(rhs.incCol),
          data(rhs.data)
    {
    }

    const ElementType &
    operator()(Index i, Index j) const
    {
        assert(i<m && j<n);
        return data[i*incRow + j*incCol];
    }

    ElementType &
    operator()(Index i, Index j)
    {
        assert(i<m && j<n);
        return data[i*incRow + j*incCol];
    }

    View
    operator()(Index i, Index j, Index numRows, Index numCols)
    {
        assert(i+numRows<=m);
        assert(j+numCols<=n);
        return View(numRows, numCols, &(operator()(i,j)), incRow, incCol);
    }

    const Index     m, n, incRow, incCol;
    ElementType*    data;
};

template <typename Matrix1, typename Matrix2>
struct GeMatrixCombinedConstView
{
    typedef typename std::common_type<typename Matrix1::ElementType,
        typename Matrix2::ElementType>::type ElementType;
    typedef typename std::common_type<typename Matrix1::Index,
        typename Matrix2::Index>::type Index;

    const Index m;
    const Index n;
    const Matrix1& A;
    const Matrix2& B;
    std::function<ElementType(ElementType, ElementType)> apply;

    template <typename ApplyOperator>
    GeMatrixCombinedConstView(const Matrix1& A, const Matrix2& B,
        ApplyOperator apply)
    : m(A.m), n(A.n), A(A), B(B), apply(apply)
    {
        assert(A.m == B.m && A.n == B.n);
    }

    const ElementType
    operator()(Index i, Index j) const
    {
        return apply(A(i, j), B(i, j));
    }

};

//
//  general template to handle elements of a matrix
//
template <typename GeMatrix, typename GeFunc>
void
applyGeMatrix(GeMatrix &A, GeFunc geFunc)
{
    for (typename GeMatrix::Index j=0; j<A.n; ++j) {
        for (typename GeMatrix::Index i=0; i<A.m; ++i) {
            geFunc(i, j, A(i,j));
        }
    }
}

//
//  Interface for bench
//
template <typename GeMatrix, typename Initvalue>
void
initGeMatrix(GeMatrix &A,Initvalue initValue)
{
    bench::initGeMatrix(A.m, A.n, A.data, A.incRow, A.incCol, initValue);
}

template <typename GeMatrix>
typename GeMatrix::ElementType
asumDiffGeMatrix(const GeMatrix &A, const GeMatrix &B)
{
    return bench::asumDiffGeMatrix(A.m, A.n,
                                   A.data, A.incRow, A.incCol,
                                   B.data, B.incRow, B.incCol);
}

} // namespace matvec

#endif // HPC_GEMATRIX_H
