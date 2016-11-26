#ifndef INC_MATRIX_CLASS_HPP
#define INC_MATRIX_CLASS_HPP

#include <cstddef> /* needed for std::size_t and std::ptrdiff_t */
#include <cassert> /* needed for assert */
#include <fmt/printf.hpp> /* needed for fmt::printf */

enum class StorageOrder {ColMajor, RowMajor};

template <typename T, typename I>
    struct GeMatrixView;

template <typename T, typename I=std::size_t>
struct GeMatrix
{
    typedef T                       ElementType;
    typedef I                       Index;
    typedef GeMatrix<T,Index>       NoView;
    typedef GeMatrixView<T,Index>   View;

    const Index m, n, incRow, incCol;
    ElementType* data;

    GeMatrix(Index m, Index n, StorageOrder order)
        :   m(m), n(n),
            incRow(order == StorageOrder::ColMajor? 1: n),
            incCol(order == StorageOrder::RowMajor? 1: m),
            data(new GeMatrix::ElementType[m*n]) {
    }

    const ElementType& operator()(Index i, Index j) const
    {
        assert(i < m && j < n);
        return data[i*incRow + j*incCol];
    }

    ElementType& operator()(Index i, Index j)
    {
        assert(i < m && j < n);
        return data[i*incRow + j*incCol];
    }

    View operator()(Index i, Index j, Index numRows, Index numCols)
    {
        assert(i+numRows<=m);
        assert(j+numCols<=n);
        return View(numRows, numCols, &(operator()(i,j)), incRow, incCol);
    }

    void init()
    {
        for (Index i = 0; i < m; ++i) {
            for (Index j = 0; j < n; ++j) {
            data[i*incRow + j*incCol] = j * n + i + 1;
            }
        }
    }

    void initRand()
    {
        for (Index j=0; j<n; ++j) {
            for (Index i=0; i<m; ++i) {
                data[i*incRow+j*incCol] = ((ElementType)rand()
                                        - RAND_MAX/2)*200/RAND_MAX;
            }
        }
    }

    void print()
    {
        for (Index i = 0; i < m; ++i) {
            fmt::printf("  ");
            for (Index j = 0; j < n; ++j) {
                fmt::printf(" %5.1lf", data[i*incRow + j*incCol]);
            }
            fmt::printf("\n");
        }
        fmt::printf("\n");
    }

    ~GeMatrix()
    {
        delete[] data;
    }
};

template <typename T, typename I=std::size_t>
struct GeMatrixView
{
    typedef T                       ElementType;
    typedef I                       Index;
    typedef GeMatrix<T,Index>       NoView;
    typedef GeMatrixView<T,Index>   View;

    const Index     m, n, incRow, incCol;
    ElementType*    data;

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
};

#endif  /* INC_MATRIX_CLASS_HPP */
