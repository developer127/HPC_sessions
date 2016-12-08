#include <thread>   //compile with -lpthread
#include <random>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/print.h>
#include <hpc/matvec/slices.h>
#include <hpc/aux/walltime.h>

template<typename MA>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
randomInit(MA& A) {
    using ElementType = typename MA::ElementType;
    using Index = typename MA::Index;

    std::random_device random;
    std::mt19937 mt(random());
    std::uniform_real_distribution<ElementType> uniform(-100,100);

    hpc::matvec::apply(A, [&](ElementType& val, Index i, Index j) -> void {
        val = uniform(mt);
    });
}

volatile double global;

int main() {
    using namespace hpc::matvec;
    using namespace hpc::aux;
    unsigned int nof_threads = std::thread::hardware_concurrency();

    WallTime<double> wall_time;

    constexpr unsigned int N = 5005;

    std::vector<std::thread> threads(nof_threads);
    wall_time.tic();
    {
        GeMatrix<double> A(N, N, StorageOrder::RowMajor);

        if(A.incRow<A.incCol)       //colmajor
        {
            Slices<GeMatrix<double>::Index> slices(nof_threads, A.numCols);

            for (unsigned int j = 0; j < nof_threads; ++j) {
                auto A_ = A(0, j*slices.offset, A.numRows, slices.size(j));
                threads[j] = std::thread( [=] () mutable {
                                randomInit(A_);
                            });
            }
        } else {                    // rowMajor

            Slices<GeMatrix<double>::Index>  slices(nof_threads, A.numRows);

            for (unsigned int i = 0; i < nof_threads; ++i) {
                auto A_ = A(i*slices.offset, 0, slices.size(i), A.numCols);
                threads[i] = std::thread( [=] () mutable {
                                randomInit(A_);
                            });
            }
        }

        for(unsigned int index = 0; index < nof_threads; ++index) {
            threads[index].join();
        }
        global = A(7, 2);
    }
    double t1 = wall_time.toc();


    wall_time.tic();
    {
        GeMatrix<double> A(N, N, StorageOrder::RowMajor);
        randomInit(A);
        global = A(7, 2);
    }
    double t2 = wall_time.toc();

    fmt::printf("t1 = %4.3lf,\nt2 = %4.3lf\n", t1, t2);

}
