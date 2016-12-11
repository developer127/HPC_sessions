#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/copy.h>
#include <hpc/matvec/print.h>
#include <hpc/matvec/randomInit.hpp>


int main() {
    using namespace hpc::matvec;
    GeMatrix<double> A(4,5);
    randomInit(A);
    GeMatrix<double> B(4,5);
    //copy(A, B);
    B = A;
    print(A);
    print(B);
    return 0;
}
