#ifndef INC_BLAS_HPP
#define INC_BLAS_HPP

#include <cassert>
#include <cstddef>
#include "matrix_class.hpp"

namespace ulmBLAS {

void
gecopy(const GeMatrix& A, GeMatrix& B);

void
geaxpy(double alpha, const GeMatrix& X, GeMatrix& Y);

void
gescal(double alpha, GeMatrix& X);

}

namespace refColMajor {

void
gemm(double alpha, const GeMatrix& A, const GeMatrix& B,
     double beta, GeMatrix& C);

}

#endif      //INC_BLAS_HPP
