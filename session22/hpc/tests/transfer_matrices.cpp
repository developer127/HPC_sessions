#include <cassert>
#include <mpi.h>
#include <printf.hpp>
#include <hpc/matvec/gematrix.h>
#include <hpc/mpi/matrix.hpp>
#include <hpc/mpi/vector.hpp>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/print.h>

int main(int argc, char** argv) {
   MPI_Init(&argc, &argv);

   int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
   int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   assert(nof_processes == 2);

   using namespace hpc::matvec;
   using namespace hpc::mpi;
   using namespace hpc::aux;

   GeMatrix<double> A(3, 7, StorageOrder::ColMajor);
   GeMatrix<double> B(3, 7, StorageOrder::RowMajor);
   using Index = GeMatrix<double>::Index;
   MPI_Datatype datatype_A = get_type(A);
   MPI_Datatype datatype_B = get_type(B);

   if (rank == 0) {
      MPI_Aint true_extent;
      MPI_Aint true_lb;
      MPI_Type_get_true_extent(datatype_A, &true_lb, &true_extent);
      MPI_Aint extent;
      MPI_Aint lb;
      MPI_Type_get_extent(datatype_A, &lb, &extent);
      std::ptrdiff_t size = sizeof(double) * A.numRows * A.numCols;
      fmt::printf("true extent of A = %d\n", true_extent);
      fmt::printf("extent of A = %d\n", extent);
      fmt::printf("size = %zd\n", size);
      assert(true_extent == size);
   }

   if (rank == 0) {
      apply(A, [](double& val, Index i, Index j) -> void {
	 val = i * 100 + j;
      });
      MPI_Send(&A(0, 0), 1, datatype_A, 1, 0, MPI_COMM_WORLD);
      MPI_Status status;
      MPI_Recv(&B(0, 0), 1, datatype_B, 1, 0, MPI_COMM_WORLD, &status);
      apply(A, [&](double& val, Index i, Index j) -> void {
	 if (val != B(i, j)) {
	    fmt::printf("verification failed for (%d,%d): %lg vs %lg\n",
	       i, j, val, B(i, j));
	 }
      });
   } else {
      MPI_Status status;
      MPI_Recv(&A(0, 0), 1, datatype_A, 0, 0, MPI_COMM_WORLD, &status);
      MPI_Send(&A(0, 0), 1, datatype_A, 0, 0, MPI_COMM_WORLD);
   }

   MPI_Finalize();
}
