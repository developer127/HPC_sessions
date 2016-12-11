#ifndef HPC_AUX_SLICE_H
#define HPC_AUX_SLICE_H 1

#include <cassert>

namespace hpc { namespace aux {

template<typename T>
struct Slices {
   Slices(T nof_threads, T problem_size) :
	 nof_threads((assert(nof_threads > 0), nof_threads)),
	 problem_size(problem_size),
	 remainder(problem_size % nof_threads),
	 slice_size(problem_size / nof_threads) {
   }
   T offset(T index) {
      assert(index < nof_threads);
      if (index < remainder) {
	 return index * (slice_size + 1);
      } else {
	 return remainder * (slice_size + 1) +
	        (index - remainder) * slice_size;
      }
   }
   T size(T index) {
      assert(index < nof_threads);
      if (index < remainder) {
	 return slice_size + 1;
      } else {
	 return slice_size;
      }
   }
   T nof_threads; T problem_size;
   T remainder; T slice_size;
};

template<typename T1, typename T2, typename Body>
void foreach_slice(T1 nof_threads, T2 problem_size, Body body) {
   using T = typename std::common_type<T1, T2>::type;
   Slices<T> slices{T(nof_threads), T(problem_size)};
   for (T index = T(0); index < T(nof_threads); ++index) {
      body(slices.offset(index), slices.size(index));
   }
}

} } // namespace aux, hpc

#endif // HPC_AUX_SLICE_H
