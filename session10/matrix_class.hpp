#include <cstddef> /* needed for std::size_t and std::ptrdiff_t */
#include <cassert> /* needed for assert */
#include <printf.hpp> /* needed for fmt::printf */

enum class StorageOrder {ColMajor, RowMajor};

template <typename T>
struct Matrix {
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

   void init() {
      for (std::size_t i = 0; i < m; ++i) {
         for (std::size_t j = 0; j < n; ++j) {
            data[i*incRow + j*incCol] = j * n + i + 1;
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
   Matrix<double> A(7, 8, StorageOrder::ColMajor);
   /*A.init();
   A(5, 3) = -1.0;
   fmt::printf("A =\n");
   A.print();
   */
}
