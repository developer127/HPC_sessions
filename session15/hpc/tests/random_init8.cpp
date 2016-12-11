#include <thread>
#include <random>
#include <vector>
#include <cstdlib>
#include <mutex>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/print.h>
#include <hpc/aux/slices.h>

struct MyRandomGenerator {
    MyRandomGenerator() : mt(std::random_device()()), uniform(-100, 100) {
    }
    double gen() {
        std::lock_guard<std::mutex> lock(mutex);
        return uniform(mt);
    }

    private:
        std::mutex mutex;
        std::mt19937 mt;
        std::uniform_real_distribution<double> uniform;
};

template<typename MA, typename RNG>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
randomInit(MA& A, RNG& myGen) {
    using ElementType = typename MA::ElementType;
    using Index = typename MA::Index;

    hpc::matvec::apply(A, [&](ElementType& val, Index i, Index j) -> void {
                       val = myGen.gen();
    });
}

int main() {
    using namespace hpc::matvec;
    using namespace hpc::aux;

    GeMatrix<double> A(51, 7);
    unsigned int nof_threads = std::thread::hardware_concurrency();

    MyRandomGenerator myGen;

    std::vector<std::thread> threads(nof_threads);
    Slices<GeMatrix<double>::Index> slices(nof_threads, A.numRows);
    for (int index = 0; index < nof_threads; ++index) {
        auto firstRow = slices.offset(index);
        auto numRows = slices.size(index);
        auto A_ = A(firstRow, 0, numRows, A.numCols);
        threads[index] = std::thread([=,&myGen]() mutable {
                                     randomInit(A_, myGen); });
    }
    for (int index = 0; index < nof_threads; ++index) {
        threads[index].join();
    }
    /* print a small block of each of the initialized matrices */
    print(A, "A");
}
