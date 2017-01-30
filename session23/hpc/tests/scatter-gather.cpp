#include <mpi.h>
#include <hpc/mpi/matrix.hpp>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/print.h>
#include <hpc/matvec/apply.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    using namespace hpc::matvec;
    using namespace hpc::mpi;

    using Matrix = GeMatrix<double>;
    using Index = typename Matrix::Index;
    int share = 3;
    int num_rows = nof_processes * share;
    int num_cols = 5;

    Matrix A(rank==0 ? num_rows : 1,
             rank==0 ? num_cols : 1);
    Matrix B(share, num_cols, StorageOrder::RowMajor); /* individual share */

    if (rank == 0) {
        apply(A, [](double& val, Index i, Index j) -> void {
            val = i * 100 + j;
        });
    }

    /* using MPI_Scatter: scatter A / receive our share into B */

    MPI_Scatter(&A(0,0), share, get_row_type(A),
                &B(0,0), share, get_row_type(B),
                0, MPI_COMM_WORLD);

    apply(B, [=](double& val, Index i, Index j) -> void {
        val += 10000 * (rank + 1);
    });

    /* using MPI_Gather: gather into A / send our share from B */

    MPI_Gather(&B(0,0), share, get_row_type(B),
               &A(0,0), share, get_row_type(A),
               0, MPI_COMM_WORLD);

    if (rank == 0) {
        print(A, "A");
    }
    MPI_Finalize();
}
