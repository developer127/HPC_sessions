#ifndef HPC_CUDA_COPY_H
#define HPC_CUDA_COPY_H 1

#include <cassert>
#include <hpc/cuda/check.h>
#include <hpc/cuda/gematrix.h>
#include <hpc/cuda/densevector.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/densevector.h>

namespace hpc { namespace cuda {

template<typename T, typename IndexX, typename IndexY>
void
copy(const hpc::matvec::DenseVector<T, IndexX>& x, DenseVector<T, IndexY>& y) {
   assert(x.length == y.length);
   CHECK_CUDA(cudaMemcpy, y.data, x.data, x.length * sizeof(T),
      cudaMemcpyHostToDevice);
}

template<typename T, typename IndexX, typename IndexY>
void
copy(const DenseVector<T, IndexX>& x, hpc::matvec::DenseVector<T, IndexY>& y) {
   assert(x.length == y.length);
   CHECK_CUDA(cudaMemcpy, y.data, x.data, x.length * sizeof(T),
      cudaMemcpyDeviceToHost);
}

template<typename T, typename IndexA, typename IndexB>
void
copy(const hpc::matvec::GeMatrix<T, IndexA>& A, GeMatrix<T, IndexB>& B) {
   assert(A.numRows == B.numRows && A.numCols == B.numCols &&
      A.incRow == B.incRow);
   CHECK_CUDA(cudaMemcpy, B.data, A.data, A.numRows * A.numCols * sizeof(T),
      cudaMemcpyHostToDevice);
}

template<typename T, typename IndexA, typename IndexB>
void
copy(const GeMatrix<T, IndexA>& A, hpc::matvec::GeMatrix<T, IndexB>& B) {
   assert(A.numRows == B.numRows && A.numCols == B.numCols &&
      A.incRow == B.incRow);
   CHECK_CUDA(cudaMemcpy, B.data, A.data, A.numRows * A.numCols * sizeof(T),
      cudaMemcpyDeviceToHost);
}

} } // namespaces cuda, hpc

#endif // HPC_CUDA_COPY_H
