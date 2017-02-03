#include <cassert>
#include <cmath>
#include <printf.hpp>
#include <mpi.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/mpi/matrix.hpp>
#include <hpc/mpi/vector.hpp>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/copy.h>
#include <hpc/matvec/matrix2pixbuf.h>
#include <hpc/aux/slices.h>
#include <hpc/aux/hsvcolor.h>
#include <unistd.h>

constexpr auto PI = std::acos(-1.0);
constexpr auto E = std::exp(1.0);
constexpr auto E_POWER_MINUS_PI = std::pow(E, -PI);

template<typename Matrix>
typename Matrix::ElementType jacobi_iteration(const Matrix& A, Matrix& B)
{
    using ElementType = typename Matrix::ElementType;
    using Index = typename Matrix::Index;
    assert(A.numRows > 2 && A.numCols > 2);
    ElementType maxdiff = 0;
    for (Index i = 1; i + 1 < B.numRows; ++i) {
        for (Index j = 1; j + 1 < B.numCols; ++j) {
            B(i, j) = 0.25
                    * (A(i - 1, j) + A(i + 1, j) + A(i, j - 1) + A(i, j + 1));
            double diff = std::fabs(A(i, j) - B(i, j));
            if (diff > maxdiff) maxdiff = diff;
        }
    }
    return maxdiff;
}

template<typename Matrix>
void exchange_with_neighbors(Matrix& A, /* ranks of the neighbors */
                             int previous, int next,
                    /* data type for an inner row, i.e. without the border */
                             MPI_Datatype rowtype)
{
    MPI_Request request[4];
    // upward
    MPI_Isend(&A(1,1), 1, rowtype,
              previous, 0, MPI_COMM_WORLD, &request[0]);
    MPI_Irecv(&A(A.numRows-1,1), 1, rowtype,
              next, 0, MPI_COMM_WORLD, &request[1]);

    // downward
    MPI_Isend(&A(A.numRows-2,1), 1, rowtype,
              next, 0, MPI_COMM_WORLD, &request[2]);
    MPI_Irecv(&A(0,1), 1, rowtype,
              previous, 0, MPI_COMM_WORLD, &request[3]);

    //wait until everything is done
    MPI_Status status;
    for(int j=0; j<4; ++j) {
        MPI_Wait(&request[j], &status);
    }
}

int main(int argc, char** argv) {
   MPI_Init(&argc, &argv);

   int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);
   int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   using namespace hpc::matvec;
   using namespace hpc::mpi;
   using namespace hpc::aux;
   using Matrix = GeMatrix<double>;

   /* initialize the entire matrix, including its borders */
   Matrix A(100, 100, StorageOrder::RowMajor);
   using Index = GeMatrix<double>::Index;
   if (rank == 0) {
      apply(A, [&](double& val, Index i, Index j) -> void {
	 if (j == 0) {
	    val = std::sin(PI * ((double)i/(A.numRows-1)));
	 } else if (j == A.numCols - 1) {
	    val = std::sin(PI * ((double)i/(A.numRows-1))) * E_POWER_MINUS_PI;
	 } else {
	    val = 0;
	 }
      });
   }

   /* we use matrices B1 and B2 to work in our set of rows */
   Slices<Index> slices(nof_processes, A.numRows - 2);
   Matrix B1(slices.size(rank) + 2, A.numCols, StorageOrder::RowMajor);
   apply(B1, [](double& val, Index i, Index j) -> void { val = 0; });
   auto B = B1(1, 0, B1.numRows - 2, B1.numCols);

   /* distribute main body of A include left and right border */
   auto A_ = A(1, 0, A.numRows - 2, A.numCols);
   scatter_by_row(A_, B, 0, MPI_COMM_WORLD);

   /* distribute first and last row of A */
   if (rank == 0) {
      copy(A(0, 0, 1, A.numCols), B1(0, 0, 1, B1.numCols));
   }
   MPI_Datatype full_rowtype = get_row_type(B);
   if (nof_processes == 1) {
      copy(A(A.numRows-1, 0, 1, A.numCols), B1(B.numRows-1, 0, 1, B1.numCols));
   } else if (rank == 0) {
      MPI_Send(&A(A.numRows-1, 0), 1, full_rowtype, nof_processes-1, 0,
	 MPI_COMM_WORLD);
   } else if (rank == nof_processes - 1) {
      MPI_Status status;
      MPI_Recv(&B1(B.numRows-1, 0), 1, full_rowtype,
	 0, 0, MPI_COMM_WORLD, &status);
   }

   Matrix B2(B1.numRows, B1.numCols, StorageOrder::RowMajor);
   copy(B1, B2); /* actually just the border needs to be copied */

   /* compute type for inner rows without the border */
   auto B_inner = B(0, 1, 1, A.numCols - 2);
   MPI_Datatype inner_rowtype = get_type(B_inner);

   int previous = rank == 0? MPI_PROC_NULL: rank-1;
   int next = rank == nof_processes-1? MPI_PROC_NULL: rank+1;

   double eps = 1e-6; unsigned int iterations;
   for (iterations = 0; ; ++iterations) {
      double maxdiff = jacobi_iteration(B1, B2);
      exchange_with_neighbors(B2, previous, next, inner_rowtype);
      maxdiff = jacobi_iteration(B2, B1);
      if (iterations % 10 == 0) {
	 double global_max;
	 MPI_Reduce(&maxdiff, &global_max, 1, get_type(maxdiff),
	    MPI_MAX, 0, MPI_COMM_WORLD);
	 MPI_Bcast(&global_max, 1, get_type(maxdiff), 0, MPI_COMM_WORLD);
	 if (global_max < eps) break;
      }
      exchange_with_neighbors(B1, previous, next, inner_rowtype);
   }
   if (rank == 0) fmt::printf("%d iterations\n", iterations);

   gather_by_row(B, A_, 0, MPI_COMM_WORLD);

   MPI_Finalize();

   if (rank == 0) {
      auto pixbuf = create_pixbuf(A, [](double val) -> HSVColor<double> {
	 return HSVColor<double>((1-val) * 240, 1, 1);
      }, 8);
      gdk_pixbuf_save(pixbuf, "jacobi.jpg", "jpeg", nullptr,
	 "quality", "100", nullptr);
   }
}
