#include <thread>   //compile with -lpthread
#include <random>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/print.h>
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
            std::size_t jobSize = (A.numCols+nof_threads-1)/ nof_threads;
            unsigned int remainder = (A.numCols%jobSize==0)? jobSize
                                                           : A.numCols%jobSize;

            for (unsigned int j = 0; j < nof_threads-1; ++j) {
                auto A_ = A(0, j*jobSize, A.numRows, jobSize);
                threads[j] = std::thread( [=] () mutable {
                                randomInit(A_);
                            });
            }
            auto A_ = A(0, (nof_threads-1)*jobSize, A.numRows, remainder);
            threads[nof_threads-1] = std::thread( [=] () mutable {
                            randomInit(A_);
                        });
        } else {                    // rowMajor

            std::size_t jobSize = (A.numRows+nof_threads-1)/ nof_threads;
            unsigned int remainder = (A.numRows%jobSize==0)? jobSize
                                                           : A.numRows%jobSize;

            for (unsigned int i = 0; i < nof_threads-1; ++i) {
                auto A_ = A(i*jobSize, 0, jobSize, A.numCols);
                threads[i] = std::thread( [=] () mutable {
                                randomInit(A_);
                            });
            }
            auto A_ = A((nof_threads-1)*jobSize, 0, remainder, A.numCols);
            threads[nof_threads-1] = std::thread( [=] () mutable {
                            randomInit(A_);
                        });

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

    /* print a small block of each of the initialized matrices
    auto A1 = A(996, 996, 5, 5);
    print(A1, "A1");
    auto A2 = A(100, 100, 5, 5);
    print(A2, "A2");
    */
}
