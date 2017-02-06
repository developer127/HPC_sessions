#ifndef HPC_MPI_MATRIX_H
#define HPC_MPI_MATRIX_H 1

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>
#include <mpi.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/isgematrix.h>
#include <hpc/matvec/copy.h>
#include <hpc/mpi/fundamental.h>
#include <hpc/aux/slices.h>

namespace hpc { namespace mpi {

template<typename Matrix>
typename std::enable_if<hpc::matvec::IsGeMatrix<Matrix>::value,
   MPI_Datatype>::type
get_row_type(const Matrix& A) {
   using ElementType = typename Matrix::ElementType;
   MPI_Datatype rowtype;
   MPI_Type_vector(
      /* count = */ A.numCols,
      /* blocklength = */ 1,
      /* stride = */ A.incCol,
      /* element type = */ get_type(A(0, 0)),
      /* newly created type = */ &rowtype);

   /* in case of row major we are finished */
   if (A.incRow == A.numCols) {
      MPI_Type_commit(&rowtype);
      return rowtype;
   }

   /* the extent of the MPI data type does not match
      the offset of subsequent rows -- this is a problem
      whenever we want to handle more than one row;
      to fix this we need to use the resize function
      which allows us to adapt the extent to A.incRow */
   MPI_Datatype resized_rowtype;
   MPI_Type_create_resized(rowtype, 0, /* lb remains unchanged */
      A.incRow * sizeof(ElementType), &resized_rowtype);
   MPI_Type_commit(&resized_rowtype);
   return resized_rowtype;
}

template<typename Matrix>
typename std::enable_if<hpc::matvec::IsGeMatrix<Matrix>::value,
   MPI_Datatype>::type
get_type(const Matrix& A) {
   using ElementType = typename Matrix::ElementType;
   MPI_Datatype datatype;
   if (A.incCol == 1) {
      MPI_Type_vector(
	 /* count = */ A.numRows,
	 /* blocklength = */ A.numCols,
	 /* stride = */ A.incRow,
	 /* element type = */ get_type(A(0, 0)),
	 /* newly created type = */ &datatype);
   } else {
      /* vector of row vectors */
      MPI_Datatype rowtype = get_row_type(A);
      MPI_Type_contiguous(A.numRows, rowtype, &datatype);
   }
   MPI_Type_commit(&datatype);
   return datatype;
}

template<typename MA, typename MB>
typename std::enable_if<
   hpc::matvec::IsGeMatrix<MA>::value && hpc::matvec::IsGeMatrix<MB>::value &&
      std::is_same<typename MA::ElementType, typename MB::ElementType>::value,
   int>::type
scatter_by_row(const MA& A, MB& B, int root, MPI_Comm comm) {
   assert(A.numCols == B.numCols);

   int nof_processes; MPI_Comm_size(comm, &nof_processes);
   int rank; MPI_Comm_rank(comm, &rank);

   hpc::aux::Slices<int> slices(nof_processes, A.numRows);
   std::vector<int> counts(nof_processes);
   std::vector<int> offsets(nof_processes);

   MPI_Datatype rowtype_A = get_row_type(A);
   for (int i = 0; i < nof_processes; ++i) {
      if (i < A.numRows) {
	 counts[i] = slices.size(i);
	 offsets[i] = slices.offset(i);
      } else {
	 counts[i] = 0; offsets[i] = 0;
      }
   }
   int recvcount = counts[rank];
   assert(B.numRows == recvcount);

   MPI_Datatype rowtype_B = get_row_type(B);

   /* OpenMPI implementation of Debian wheezy expects void* instead
      of const void*; hence we need to remove const */
   return MPI_Scatterv((void*) &A(0, 0), &counts[0], &offsets[0], rowtype_A,
      &B(0, 0), recvcount, rowtype_B, root, comm);
}

template<typename MA, typename MB>
typename std::enable_if<
   hpc::matvec::IsGeMatrix<MA>::value && hpc::matvec::IsGeMatrix<MB>::value &&
      std::is_same<typename MA::ElementType, typename MB::ElementType>::value,
   int>::type
gather_by_row(const MA& A, MB& B, int root, MPI_Comm comm) {
   assert(A.numCols == B.numCols);

   int nof_processes; MPI_Comm_size(comm, &nof_processes);
   int rank; MPI_Comm_rank(comm, &rank);

   hpc::aux::Slices<int> slices(nof_processes, B.numRows);
   std::vector<int> counts(nof_processes);
   std::vector<int> offsets(nof_processes);

   for (int i = 0; i < nof_processes; ++i) {
      if (i < B.numRows) {
	 counts[i] = slices.size(i);
	 offsets[i] = slices.offset(i);
      } else {
	 counts[i] = 0; offsets[i] = 0;
      }
   }
   int sendcount = counts[rank];
   assert(A.numRows == sendcount);

   MPI_Datatype rowtype_A = get_row_type(A);
   MPI_Datatype rowtype_B = get_row_type(B);

   /* OpenMPI implementation of Debian wheezy expects void* instead
      of const void*; hence we need to remove const */
   return MPI_Gatherv((void*) &A(0, 0), sendcount, rowtype_A,
      &B(0, 0), &counts[0], &offsets[0], rowtype_B, root, comm);
}

template<typename MA, typename MB>
typename std::enable_if<
   hpc::matvec::IsGeMatrix<MA>::value && hpc::matvec::IsGeMatrix<MB>::value &&
      std::is_same<typename MA::ElementType, typename MB::ElementType>::value,
   int>::type
scatter_by_block(const MA& A, MB& B, int root,
      MPI_Comm grid, int dims[2], int coords[2], int overlap = 0) {
   assert(overlap < A.numRows && overlap < A.numCols);

   int nof_processes; MPI_Comm_size(grid, &nof_processes);
   int rank; MPI_Comm_rank(grid, &rank);
   hpc::aux::Slices<int> rows(dims[0], A.numRows - 2*overlap);
   hpc::aux::Slices<int> columns(dims[1], A.numCols - 2*overlap);

   if (rank == root) {
      MPI_Request requests[nof_processes-1]; int ri = 0;
      for (int i = 0; i < nof_processes; ++i) {
	 int coords_[2];
	 MPI_Cart_coords(grid, i, 2, coords_);
	 auto A_ = A(rows.offset(coords_[0]), columns.offset(coords_[1]),
	    rows.size(coords_[0]) + 2*overlap,
	    columns.size(coords_[1]) + 2*overlap);
	 if (i == root) {
	    hpc::matvec::copy(A_, B);
	 } else {
	    MPI_Isend(&A_(0, 0), 1, get_type(A_),
	       i, 0, grid, &requests[ri++]);
	 }
      }
      for (auto& request: requests) {
	 MPI_Status status;
	 MPI_Wait(&request, &status);
      }
   } else {
      MPI_Status status;
      MPI_Recv(&B(0, 0), 1, get_type(B), root, 0, grid, &status);
   }
}

template<typename MA, typename MB>
typename std::enable_if<
   hpc::matvec::IsGeMatrix<MA>::value && hpc::matvec::IsGeMatrix<MB>::value &&
      std::is_same<typename MA::ElementType, typename MB::ElementType>::value,
   int>::type
gather_by_block(const MA& A, MB& B, int root,
      MPI_Comm grid, int dims[2], int coords[2], int overlap = 0) {
   assert(overlap < A.numRows && overlap < A.numCols);

   int nof_processes; MPI_Comm_size(grid, &nof_processes);
   int rank; MPI_Comm_rank(grid, &rank);
   hpc::aux::Slices<int> rows(dims[0], B.numRows - 2*overlap);
   hpc::aux::Slices<int> columns(dims[1], B.numCols - 2*overlap);

   auto A_ = A(overlap, overlap, A.numRows - 2*overlap, A.numCols - 2*overlap);
   if (rank == root) {
      MPI_Request requests[nof_processes-1]; int ri = 0;
      for (int i = 0; i < nof_processes; ++i) {
	 int coords_[2];
	 MPI_Cart_coords(grid, i, 2, coords_);
	 auto B_ = B(rows.offset(coords_[0]) + overlap,
	    columns.offset(coords_[1]) + overlap,
	    rows.size(coords_[0]), columns.size(coords_[1]));
	 if (i == root) {
	    hpc::matvec::copy(A_, B_);
	 } else {
	    MPI_Irecv(&B_(0, 0), 1, get_type(B_),
	       i, 0, grid, &requests[ri++]);
	 }
      }
      for (auto& request: requests) {
	 MPI_Status status;
	 MPI_Wait(&request, &status);
      }
   } else {
      MPI_Send(&A_(0, 0), 1, get_type(A_), root, 0, grid);
   }
}

} } // namespaces mpi, hpc

#endif // HPC_MPI_MATRIX_H
