#ifndef HPC_AUX_REPOOL_H
#define HPC_AUX_REPOOL_H 1

#include <cassert>
#include <condition_variable>
#include <cstdlib>
#include <mutex>
#include <random>
#include <vector>

namespace hpc { namespace aux {

template<typename T>
struct RandomEnginePool {
      using EngineType = T;
      RandomEnginePool(std::size_t size) :
	    size(size), nof_free_engines(size),
	    inuse(size), engines(size) {
	 std::random_device r;
	 for (std::size_t i = 0; i < size; ++i) {
	    engines[i].seed(r()); inuse[i] = false;
	 }
      }
      T& get() {
	 std::unique_lock<std::mutex> lock(mutex);
	 if (nof_free_engines == 0) {
	    cv.wait(lock);
	 }
	 for (std::size_t i = 0; i < size; ++i) {
	    if (!inuse[i]) {
	       inuse[i] = true; --nof_free_engines;
	       return engines[i];
	    }
	 }
	 assert(false);
      }
      void free(T& engine) {
	 {
	    std::unique_lock<std::mutex> lock(mutex);
	    bool found = false;
	    for (std::size_t i = 0; i < size; ++i) {
	       if (&engine == &engines[i]) {
		  inuse[i] = false; ++nof_free_engines;
		  found = true; break;
	       }
	    }
	    assert(found);
	 }
	 cv.notify_one();
      }
   private:
      std::mutex mutex;
      std::condition_variable cv;
      std::size_t size;
      std::size_t nof_free_engines;
      std::vector<bool> inuse;
      std::vector<T> engines;
};

template<typename T>
struct RandomEngineGuard {
   using EngineType = T;
   RandomEngineGuard(RandomEnginePool<T>& pool) :
      pool(pool), engine(pool.get()) {
   }
   ~RandomEngineGuard() {
      pool.free(engine);
   }
   T& get() {
      return engine;
   }
   RandomEnginePool<T>& pool;
   T& engine;
};

} } // namespace aux, hpc

#endif // HPC_AUX_REPOOL_H
