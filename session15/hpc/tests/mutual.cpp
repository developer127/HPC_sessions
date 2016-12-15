#include <mutex>
#include <thread>
#include <vector>
#include <fmt/printf.hpp>
#include <hpc/aux/walltime.h>

class Counter {
   public:
      using Value = unsigned long long int;

      Counter() : counter(0) {
      }
      unsigned long long int increment() {
	 std::lock_guard<std::mutex> lock(mutex);
	 return ++counter;
      }
   private:
      std::mutex mutex;
      Value counter;
};

int main() {
   constexpr unsigned int nof_threads = 2;
   constexpr Counter::Value max_counter = 1LL << 24;

   std::vector<std::thread> threads(nof_threads);
   hpc::aux::WallTime<double> wall_time;
   Counter counter;
   wall_time.tic();
      for (auto& t: threads) {
	 t = std::thread([&]{
	    while (counter.increment() < max_counter)
	       ;
	 });
      }
      for (auto& t: threads) {
	 t.join();
      }
   auto t = wall_time.toc();
   fmt::printf("avg time per lock = %.2lf ns\n",
      t / max_counter * 1000000000LL);
}
