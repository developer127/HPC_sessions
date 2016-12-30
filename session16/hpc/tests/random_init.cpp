#include <thread>
#include <random>
#include <vector>
#include <cstdlib>
#include <mutex>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/print.h>
#include <hpc/aux/slices.h>
#include <hpc/aux/randomEnginePool.hpp>
#include <hpc/aux/randomEngineGuard.hpp>
#include <hpc/aux/tpool.hpp>
#include <hpc/aux/jobstatus.hpp>


template<typename MA, typename REPool, typename ThreadPool>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
randomInit(MA& A, REPool& pool, ThreadPool& tpool) {
    using ElementType = typename MA::ElementType;
    using Index = typename MA::Index;
    using EngineType = typename REPool::EngineType;

    std::uniform_real_distribution<double> uniform(-100,100);

    hpc::aux::RandomEngineGuard<EngineType> guard(pool);
    auto& engine(guard.get());
    unsigned int nof_threads = tpool.get_nof_threads();     //open race?

    hpc::aux::Slices<Index> slices(nof_threads, A.numRows);
    std::vector<hpc::aux::JobStatus> jobStatus(nof_threads);
    for (Index index = 0; index < nof_threads; ++index) {
        auto firstRow = slices.offset(index);
        auto numRows = slices.size(index);
        auto A_ = A(firstRow, 0, numRows, A.numCols);
        tpool.submit([&]() mutable {
            hpc::matvec::apply(A_, [&](ElementType& val, Index i, Index j)
                               -> void {
                       val = uniform(engine);
        });});
        jobStatus[index].is_finished();
    }

    for (auto &js: jobStatus) {
        js.wait();
    }
}

int main() {
    using namespace hpc::matvec;
    using namespace hpc::aux;

    unsigned int nof_threads = std::thread::hardware_concurrency();

    ThreadPool tpool(nof_threads);
    RandomEnginePool<std::mt19937> repool;

    GeMatrix<double> A(51, 7);

    randomInit(A, repool, tpool);

    /* print a small block of each of the initialized matrices */
    print(A, "A");
}
