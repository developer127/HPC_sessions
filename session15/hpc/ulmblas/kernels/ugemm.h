#ifndef HPC_ULMBLAS_KERNELS_UGEMM_H
#define HPC_ULMBLAS_KERNELS_UGEMM_H 1

#include <hpc/ulmblas/config.h>
#include <hpc/ulmblas/kernels/ugemm_ref.h>
#include <hpc/ulmblas/kernels/ugemm_buf.h>

#ifdef SSE
#include <hpc/ulmblas/kernels/ugemm_sse4x4.h>
#elif  AVX
#include <hpc/ulmblas/kernels/ugemm_avx4x8.h>
#elif  FMA
#include <hpc/ulmblas/kernels/ugemm_fma4x12.h>
#endif

#endif // HPC_ULMBLAS_KERNELS_UGEMM_H
