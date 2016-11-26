#ifndef INC_TEST_FUNC_HPP
#define INC_TEST_FUNC_HPP

#include <cassert>
#include <cmath>
#include <cstddef>
#include "matrix_class.hpp"
#include <fmt/printf.hpp> /* needed for fmt::printf */


namespace test{

double asumDiffGeMatrix(const GeMatrix& A, const GeMatrix& B);

void printGeMatrixInMemory(GeMatrix& A);

} // test
#endif
