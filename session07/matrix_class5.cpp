#include <cstddef> /* needed for std::size_t and std::ptrdiff_t */
#include <cassert> /* needed for assert */
//#include <printf.hpp> /* needed for fmt::printf */
#include <fmt/format.h> /* needed for fmt::printf */

enum class StorageOrder {ColMajor, RowMajor};

struct Matrix {
   const std::size_t m; /* number of rows */
   const std::size_t n; /* number of columns */
   const std::ptrdiff_t incRow;
   const std::ptrdiff_t incCol;
   double* data;

   Matrix(std::size_t m, std::size_t n, StorageOrder order) :
         m(m), n(n),
         incRow(order == StorageOrder::ColMajor? 1: n),
         incCol(order == StorageOrder::RowMajor? 1: m),
         data(new double[m*n]) {
   }

   const double& operator()(std::size_t i, std::size_t j) const {
      assert(i < m && j < n);
      return data[i*incRow + j*incCol];
   }

   double& operator()(std::size_t i, std::size_t j) {
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

   ~Matrix(){
       fmt::printf("Destructor used!\n");
       delete[] data;
   }
};

void
copyMatrix(const Matrix& A, Matrix& B)
{
    assert(B.n<=A.n||B.m<=A.m);
    if(B.incRow < B.incCol) {
        for(std::size_t j = 0; j<B.n; ++j) {
            for(std::size_t i =0; i<B.m; ++i) {
                B(i,j) = A(i,j);
            }
        }
    } else {
        for(std::size_t i =0; i<B.m; ++i) {
            for(std::size_t j = 0; j<B.n; ++j) {
                B(i,j) = A(i,j);
            }
        }
    }
}

int main() {
   Matrix A(7, 8, StorageOrder::ColMajor);
   Matrix B(5, 5, StorageOrder::ColMajor);
   A.init();
   copyMatrix(A,B);
   A(9, 3) = -1;
   fmt::printf("A =\n");
   A.print();
   fmt::printf("B =\n");
   B.print();
}
