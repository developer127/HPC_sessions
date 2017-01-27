#include <cassert>
#include <cstdlib>
#include <mpi.h>
#include <printf.hpp>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/densevector.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/mpi/sendvec.hpp>
#include <hpc/mpi/receivevec.hpp>

int main(int argc, char** argv) {
   MPI_Init(&argc, &argv);

   int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
   int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   assert(nof_processes == 2);

   using namespace hpc::matvec;
   using namespace hpc::mpi;
   using namespace hpc::aux;

   std::size_t nof_rows = 3;
   std::size_t nof_cols = 7;

   if (rank == 0) {
      GeMatrix<double> A(nof_rows, nof_cols, StorageOrder::RowMajor);
      using Index = GeMatrix<double>::Index;
      apply(A, [](double& val, Index i, Index j) -> void {
	 val = i * 100 + j;
      });
      auto row = A.row(2);
      auto col = A.col(0);

      sendVec(row, 1, 0);
      sendVec(col, 1, 0);

      /* receive it back for verification */
      DenseVector<double> vec1(nof_cols), vec2(nof_rows);
      receiveVec(vec1, 1, 0);
      receiveVec(vec2, 1, 0);

      /* verify it */
      apply(vec1, [=](double& val, Index i) {
	 if (val != row(i)) {
	    fmt::printf("verification failed for row(%d): %lg vs %lg\n",
	       (int)i, val, row(i));
	 }
      });
      apply(vec2, [=](double& val, Index i) {
	 if (val != col(i)) {
	    fmt::printf("verification failed for col(%d): %lg vs %lg\n",
	       (int)i, val, col(i));
	 }
      });
   } else {
      DenseVector<double> vec1(nof_cols), vec2(nof_rows);
      receiveVec(vec1, 0, 0);
      receiveVec(vec2, 0, 0);

      /* send it back for verification */

      sendVec(vec1, 0, 0);
      sendVec(vec2, 0, 0);

   }
   MPI_Finalize();
}
