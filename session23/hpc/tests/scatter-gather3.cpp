#include <mpi.h>
#include <hpc/mpi/matrix.hpp>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/print.h>
#include <hpc/matvec/apply.h>
#include <hpc/aux/slices.h>

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int nof_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    using namespace hpc::matvec;
    using namespace hpc::mpi;
    using namespace hpc::aux;

    using Matrix = GeMatrix<double>;
    using Index = typename Matrix::Index;
    int share = 3;
    int num_rows = nof_processes * share + 1;
    int num_cols = 5;

    Matrix A(num_rows, num_cols); /* entire matrix */

    Slices<int> slices(nof_processes, num_rows);
    /* individual share */
    Matrix B(slices.size(rank), num_cols, StorageOrder::RowMajor);

    if (rank == 0) {
        apply(A, [](double& val, Index i, Index j) -> void {
            val = i * 100 + j;
        });
    }

    scatter_by_row(A, B, 0, MPI_COMM_WORLD);
    apply(B, [=](double& val, Index i, Index j) -> void {
        val += 10000 * (rank + 1);
    });
    gather_by_row(B, A, 0, MPI_COMM_WORLD);

    MPI_Finalize();

    if (rank == 0) {
        print(A, "A");
    }
}
