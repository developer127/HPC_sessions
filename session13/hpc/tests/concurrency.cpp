#include <thread>
#include <fmt/printf.hpp>

int main() {
   fmt::printf("hardware concurrency = %u\n",
      std::thread::hardware_concurrency());
}
