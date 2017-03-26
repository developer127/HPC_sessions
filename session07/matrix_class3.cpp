#include <cstddef> /* needed for std::size_t and std::ptrdiff_t */
#include <printf.hpp> /* needed for fmt::printf */

enum class StorageOrder {ColMajor, RowMajor};

struct Matrix {
    const std::size_t m; /* number of rows */
    const std::size_t n; /* number of columns */
    const std::ptrdiff_t incRow;
    const std::ptrdiff_t incCol;
    double* data;

    Matrix(std::size_t m, std::size_t n, StorageOrder order) :
            m(m), n(n),
            incRow((order == StorageOrder::ColMajor)? 1 : n),
            incCol((order == StorageOrder::ColMajor)? m : 1),
            data(new double[m*n]) {
    }

    void init() {
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                data[i*incRow + j*incCol] = double(j * m + i + 1);
            }
         }
    }

    void print() {
        for (std::size_t i = 0; i < m; ++i) {
            fmt::printf("  ");
            for (std::size_t j = 0; j < n; ++j) {
                fmt::printf(" %4.1lf", data[i*incRow + j*incCol]);
            }
            fmt::printf("\n");
        }
    }
};

int main() {
    Matrix A(7, 8, StorageOrder::RowMajor);
    A.init();
    fmt::printf("A =\n"); A.print();
}
