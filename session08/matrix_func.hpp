#ifndef INC_MATRIX_FUNC_HPP
#define INC_MATRIX_FUNC_HPP 1

#include <cstddef> /* needed for std::size_t and std::ptrdiff_t */
#include <cstdlib>

void
initMatrix(std::size_t m, std::size_t n,
           double *A, std::size_t incRowA, std::size_t incColA);

#endif  /* INC_MATRIX_FUNC_HPP */
