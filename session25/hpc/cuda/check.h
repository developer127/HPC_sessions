#ifndef HPC_CUDA_CHECK_H
#define HPC_CUDA_CHECK_H 1

#ifndef __CUDACC__
#	error This source must be compiled using nvcc
#endif

/* this header file defines a macro that can be used
   to wrap all CUDA library calls. Example:

      int device; cudaGetDevice(&device);

   can be replaced by

      CHECK_CUDA(cudaGetDevice, &device);

   Whenever an error is detected, CHECK_CUDA
   aborts the execution after having printed
   a helpful message pointing to the failed
   CUDA library call.
*/

#include <cstdio>
#include <cstdlib>

namespace hpc { namespace cuda {

inline void check_cuda_error(const char* cudaop, const char* source,
       unsigned int line, cudaError_t error) {
   if (error != cudaSuccess) {
      fprintf(stderr, "%s at %s:%u failed: %s\n",
	 cudaop, source, line, cudaGetErrorString(error));
      exit(1);
   }
}

} } // namespaces cuda, hpc

#define CHECK_CUDA(opname, ...) \
   hpc::cuda::check_cuda_error(#opname, __FILE__, __LINE__, opname(__VA_ARGS__))

#endif
