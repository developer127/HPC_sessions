#include <thread>
#include <fmt/printf.hpp>
#include <hpc/aux/walltime.h>

int main() {
   unsigned int counter = 0;
   unsigned int nof_threads = 1<<15;

   hpc::aux::WallTime<double> wall_time;
   wall_time.tic();
      for (unsigned int i = 0; i < nof_threads; ++i) {
         auto t = std::thread([&]() { ++counter; });
         t.join();
      }
   auto t = wall_time.toc();
   fmt::printf("avg time per thread creation = %.2f us\n",
      t / nof_threads * 1000000L);
}
