#ifndef INC_BLAS_HPP
#define INC_BLAS_HPP

#include <cassert>
#include <cstddef>

namespace ulmBLAS {

void
gecopy(std::size_t m, std::size_t n,
       const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
       double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB);

void
geaxpy(std::size_t m, std::size_t n,
       double alpha,
       const double *X, std::ptrdiff_t incRowX, std::ptrdiff_t incColX,
       double       *Y, std::ptrdiff_t incRowY, std::ptrdiff_t incColY);

void
gescal(std::size_t m, std::size_t n,
       double alpha,
       double *X, std::ptrdiff_t incRowX, std::ptrdiff_t incColX);

void
ger(std::size_t m, std::size_t n,
    double alpha,
    const double *x, std::ptrdiff_t incX,
    const double *y, std::ptrdiff_t incY,
    double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC);


} // namespace ulmBLAS

namespace refColMajor {

void
gemm(std::size_t m, std::size_t n, std::size_t k,
     double alpha,
     const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
     const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
     double beta,
     double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC);


} // namespace refColMajor

namespace simpleBuffer {

void
gemm(std::size_t m, std::size_t n, std::size_t k,
     double alpha,
     const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
     const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
     double beta,
     double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC);

} // namespace simpleBuffer

namespace blocked {

void
pack_A(std::size_t mc, std::size_t kc,
       const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
       double *p);

void
pack_B(std::size_t kc, std::size_t nc,
       const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
       double *p);

void
ugemm(std::size_t kc, double alpha,
      const double *A, const double *B,
      double beta,
      double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC);

void
mgemm(std::size_t mc, std::size_t nc, std::size_t kc,
      double alpha,
      const double *A, const double *B,
      double beta,
      double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC);

void
gemm(std::size_t m, std::size_t n, std::size_t k,
     double alpha,
     const double *A, std::ptrdiff_t incRowA, std::ptrdiff_t incColA,
     const double *B, std::ptrdiff_t incRowB, std::ptrdiff_t incColB,
     double beta,
     double *C, std::ptrdiff_t incRowC, std::ptrdiff_t incColC);

} // namespace blocked

#endif      //INC_BLAS_HPP
