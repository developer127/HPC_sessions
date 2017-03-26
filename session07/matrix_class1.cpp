#include <cstddef> /* needed for std::size_t and std::ptrdiff_t */
#include <printf.hpp> /* needed for fmt::printf */


struct Matrix {
    std::size_t m; /* number of rows */
    std::size_t n; /* number of columns */
    std::ptrdiff_t incRow;
    std::ptrdiff_t incCol;
    double* data;

    void init()
    {
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                data[i*incRow + j*incCol] = j * m + i + 1;
            }
        }
    }
    void print_matrix() const
    {
        for (std::size_t i = 0; i < m; ++i) {
            fmt::printf("  ");
            for (std::size_t j = 0; j < n; ++j) {
                fmt::printf(" %4.1lf", data[i*incRow + j*incCol]);
            }
            fmt::printf("\n");
        }
    }
};


void transpose(Matrix& A)
{
    std::swap(A.m, A.n);
    std::swap(A.incRow, A.incCol);
}



int main() {
   Matrix A;
   A.m = 7; A.n = 8;
   A.data = new double[A.m * A.n];
   A.incRow = 1; A.incCol = 7;
   A.init();
   fmt::printf("A =\n");
   A.print_matrix();
   transpose(A);
   fmt::printf("A^t =\n");
   A.print_matrix();
}
