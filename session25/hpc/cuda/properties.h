#ifndef HPC_CUDA_PROPERTIES_H
#define HPC_CUDA_PROPERTIES_H 1

#ifndef __CUDACC__
#	error This source must be compiled using nvcc
#endif

#include <cstdio>
#include <cstdlib>
#include <hpc/cuda/check.h>

namespace hpc { namespace cuda {

inline const cudaDeviceProp& get_properties(int device = -1) {
   if (device < 0) {
      CHECK_CUDA(cudaGetDevice, &device);
   }
   static int current_device = -1;
   static cudaDeviceProp properties;
   if (device != current_device) {
      CHECK_CUDA(cudaGetDeviceProperties, &properties, device);
      current_device = device;
   }
   return properties;
}

int get_max_threads_per_block(int device = -1) {
   const cudaDeviceProp& properties(get_properties(device));
   return properties.maxThreadsPerBlock;
}

int get_number_of_multiprocessors(int device = -1) {
   const cudaDeviceProp& properties(get_properties(device));
   return properties.multiProcessorCount;
}

int get_warp_size(int device = -1) {
   const cudaDeviceProp& properties(get_properties(device));
   return properties.warpSize;
}

} } // namespaces cuda, hpc

#endif
