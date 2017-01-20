#include <mpi.h>
#include <printf.hpp>

int main(int argc, char** argv) {
   MPI_Init(&argc, &argv);

   int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

   if (rank) {
      MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
   } else {
      for (int i = 0; i + 1 < nof_processes; ++i) {
	 MPI_Status status;
	 int msg;
	 MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE,
	    0, MPI_COMM_WORLD, &status);
	 int count;
	 MPI_Get_count(&status, MPI_INT, &count);
	 if (count == 1) {
	    fmt::printf("%d\n", msg);
	 }
      }
   }

   MPI_Finalize();
}
