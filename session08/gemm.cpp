#include "benchmark.hpp"
#include "test_func.hpp"
#include "matrix_class.hpp"
#include "blas.hpp"

#ifndef COLMAJOR
#define COLMAJOR 1
#endif

#ifndef MAXDIM_M
#define MAXDIM_M    4000
#endif

#ifndef MAXDIM_N
#define MAXDIM_N    4000
#endif

#ifndef MAXDIM_K
#define MAXDIM_K    4000
#endif

#ifndef MIN_M
#define MIN_M   100
#endif

#ifndef MIN_N
#define MIN_N   100
#endif

#ifndef MIN_K
#define MIN_K   100
#endif

#ifndef MAX_M
#define MAX_M   7000
#endif

#ifndef MAX_N
#define MAX_N   7000
#endif

#ifndef MAX_K
#define MAX_K   7000
#endif

#ifndef INC_M
#define INC_M   100
#endif

#ifndef INC_N
#define INC_N   100
#endif

#ifndef INC_K
#define INC_K   100
#endif

#ifndef ALPHA
#define ALPHA   1.5
#endif

#ifndef BETA
#define BETA    1.5
#endif


int main()
{
    StorageOrder storage = StorageOrder::RowMajor;

    if(COLMAJOR==1) {
        storage = StorageOrder::ColMajor;
    }
    using namespace ulmBLAS;
    using namespace bench;
    using namespace test;

    //GeMatrix A(MAXDIM_M*MAXDIM_K);
    GeMatrix A(4,4,storage);
    GeMatrix B(4,4,storage);
    GeMatrix C(4,4,storage);
    //GeMatrix B(MAXDIM_K,MAXDIM_N,storage);

    GeMatrix C1(MAXDIM_M,MAXDIM_N,storage);
    GeMatrix C2(MAXDIM_M,MAXDIM_N,storage);
    GeMatrix C3(MAXDIM_M,MAXDIM_N,storage);

    A.init();
    A.print();
    gecopy(A,B);
    B.print();
    gescal(0,C);
    refColMajor::gemm(1,A,B,1,C);
    C.print();
    fmt::printf("diffA,B: %lf\n", asumDiffGeMatrix(A,B));
    fmt::printf("diffA,C: %lf\n", asumDiffGeMatrix(A,C));

    fmt::printf("A in Memory\n");
    printGeMatrixInMemory(A);


    return 0;
}
