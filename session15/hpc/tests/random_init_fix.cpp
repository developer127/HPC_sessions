#include <thread>
#include <random>
#include <vector>
#include <cstdlib>
#include <mutex>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/print.h>
#include <hpc/aux/slices.h>
#include <hpc/aux/randomEnginePoolFixSize.hpp>
#include <hpc/aux/randomEngineGuardFixSize.hpp>


template<typename MA, typename REPool>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
randomInit(MA& A, REPool& pool) {
    using ElementType = typename MA::ElementType;
    using Index = typename MA::Index;
    using EngineType = typename REPool::EngineType;

    std::uniform_real_distribution<double> uniform(-100,100);

    hpc::aux::RandomEngineGuard<EngineType> guard(pool);
    auto& engine(guard.get());

    hpc::matvec::apply(A, [&](ElementType& val, Index i, Index j) -> void {
                       val = uniform(engine);
    });
}

int main() {
    using namespace hpc::matvec;
    using namespace hpc::aux;

    GeMatrix<double> A(51, 7);
    unsigned int nof_threads = std::thread::hardware_concurrency();

    RandomEnginePool<std::mt19937> repool(2);

    typedef GeMatrix<double>::Index Index;

    std::vector<std::thread> threads(nof_threads);
    Slices<Index> slices(nof_threads, A.numRows);
    for (Index index = 0; index < nof_threads; ++index) {
        auto firstRow = slices.offset(index);
        auto numRows = slices.size(index);
        auto A_ = A(firstRow, 0, numRows, A.numCols);
        threads[index] = std::thread([=,&repool]() mutable {
                                     randomInit(A_, repool); });
    }
    for (Index index = 0; index < nof_threads; ++index) {
        threads[index].join();
    }
    /* print a small block of each of the initialized matrices */
    print(A, "A");
}
