#include <cstddef> /* needed for std::size_t and std::ptrdiff_t */
#include <cassert> /* needed for assert */
#include <gmpxx.h>
#include <fmt/printf.hpp> /* needed for fmt::printf */

enum class StorageOrder {ColMajor, RowMajor};

template <typename T>
struct Matrix {
    typedef T Element;

    const std::size_t m; /* number of rows */
    const std::size_t n; /* number of columns */
    const std::ptrdiff_t incRow;
    const std::ptrdiff_t incCol;
    T* data;

    Matrix(std::size_t m, std::size_t n, StorageOrder order) :
        m(m), n(n),
        incRow(order == StorageOrder::ColMajor? 1: n),
        incCol(order == StorageOrder::RowMajor? 1: m),
        data(new T[m*n]) {
    }

    const T& operator()(std::size_t i, std::size_t j) const {
        assert(i < m && j < n);
        return data[i*incRow + j*incCol];
    }

    T& operator()(std::size_t i, std::size_t j) {
        assert(i < m && j < n);
        return data[i*incRow + j*incCol];
    }

    Matrix(const Matrix& other) = delete;
    Matrix& operator=(const Matrix& other) = delete;

    ~Matrix() {
        delete[] data;
    }
};

template<typename T>
struct MatrixView {
    typedef T Element;

    const std::size_t m; /* number of rows */
    const std::size_t n; /* number of columns */
    const std::ptrdiff_t incRow;
    const std::ptrdiff_t incCol;
    T* data;

    MatrixView(std::size_t m, std::size_t n,
               T* data,
               std::ptrdiff_t incRow, std::ptrdiff_t incCol) :
        m(m), n(n),
        incRow(incRow), incCol(incCol),
        data(data){
    }

    const T& operator()(std::size_t i, std::size_t j) const {
        assert(i < m && j < n);
        return data[i*incRow + j*incCol];
    }

    T& operator()(std::size_t i, std::size_t j) {
        assert(i < m && j < n);
        return data[i*incRow + j*incCol];
    }
};



template<typename Matrix>
void init_Matrix(Matrix &A) {
    for (std::size_t i = 0; i < A.m; ++i) {
        for (std::size_t j = 0; j < A.n; ++j) {
            A(i,j) = j * A.m + i + 1;
        }
    }
}

/** implementation of print_value for diefferent data types
  * with a fallback to long double
  */
void print_value(long double value) {
    std::printf(" %4.1Lf", value);
}

void print_value(double value) {
    std::printf(" %4.1lf", value);
}

void print_value(float value) {
    std::printf(" %4.1f", value);
}

void print_value(mpq_class value) {
    std::printf(" %9s", value.get_str().c_str());
}

template<typename T>        //fallback
void print_value(T& value)
{
    std::printf(" %4.1LF", (long double) value);
}

template<typename Matrix>
void print_Matrix(const Matrix &A)
{
    for (std::size_t i = 0; i < A.m; ++i) {
        printf("  ");
        for (std::size_t j = 0; j < A.n; ++j) {
            print_value(A(i,j));
        }
        printf("\n");
    }
}

template<typename Matrix>
MatrixView<typename Matrix::Element>create_view(Matrix &A,
                   const std::size_t i, const std::size_t j,
                   const std::size_t numRows, const std::size_t numCols)
{
    assert(i+numRows <= A.m && j+numCols <= A.n);
    return MatrixView<typename Matrix::Element>(numRows, numCols,
                                                &A(i,j),
                                                A.incRow, A.incCol);
}

template<typename Matrix>
void scale_Matrix(Matrix &A, typename Matrix::Element beta)
{
    for(std::size_t i=0; i<A.m; ++i) {
        for(std::size_t j=0; j< A.n; ++j) {
            A(i,j) *= beta;
        }
    }
}

int main() {
    Matrix<mpq_class> A(7, 8, StorageOrder::ColMajor);
    init_Matrix(A);
    mpq_class a(1,5);
    scale_Matrix(A, a);
    A(5, 3) = -1;
    MatrixView<mpq_class> B(create_view(A,2,2,3,3));
    B(1,1) = 0;
    MatrixView<mpq_class> C(create_view(B,1,1,2,2));
    C(1,1) = 0;
    fmt::printf("A =\n");
    print_Matrix(A);
    fmt::printf("B =\n");
    print_Matrix(B);
    fmt::printf("C =\n");
    print_Matrix(C);
    init_Matrix(C);
    fmt::printf("C =\n");
    print_Matrix(C);
    fmt::printf("A =\n");
    print_Matrix(A);

}
