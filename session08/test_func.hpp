#ifndef INC_TEST_FUNC_HPP
#define INC_TEST_FUNC_HPP

#include <cassert>
#include <cmath>
#include <cstddef>
#include <printf.hpp> /* needed for fmt::printf */


namespace test{

double
asumDiffMatrix(std::size_t m, std::size_t n,
               const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
               double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB);

void
printGeMatrixInMemory(std::size_t m, std::size_t n,
                      const double *A);

} // test
#endif
