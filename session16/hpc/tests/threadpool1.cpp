#include <hpc/aux/tpool.hpp>
#include <fmt/printf.hpp>

int main() {
    using namespace hpc::aux;

    unsigned int nof_threads = std::thread::hardware_concurrency();
    fmt::printf("We use %d treads.\n", nof_threads);

    ThreadPool pool(nof_threads);
    /*
    pool.submit([]() { fmt::printf("Hi, this is your first job!\n"); });
    pool.submit([]() { fmt::printf(std::cout, "Now you got another job.\n"); });
    */
    pool.submit([&]() {
        fmt::printf("Hi, this is thread %x\n", std::this_thread::get_id());
        pool.submit([&]() {
            pool.submit([]() {
                fmt::printf("Hello guys, I'm %x\n", std::this_thread::get_id());
            });
            fmt::printf("Hi guys, I'm %x\n", std::this_thread::get_id());
        });
        pool.submit([&]() {
            pool.submit([]() {
                fmt::printf("O wonder, I'm %x\n", std::this_thread::get_id());
            });
            fmt::printf("Huhu guys, I'm %x\n", std::this_thread::get_id());
        });
        fmt::printf("Hi, this is again thread %x\n",
                    std::this_thread::get_id());
    });
}
