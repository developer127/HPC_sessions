#ifndef HPC_MPI_MATRIX_H
#define HPC_MPI_MATRIX_H 1

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>
#include <mpi.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/isgematrix.h>
#include <hpc/mpi/fundamental.h>
#include <hpc/aux/slices.h>

namespace hpc { namespace mpi {

template<typename Matrix>
typename std::enable_if<hpc::matvec::IsGeMatrix<Matrix>::value,
   MPI_Datatype>::type
get_row_type(const Matrix& A)
{
    using ElementType = typename Matrix::ElementType;
    MPI_Datatype rowtype;
    MPI_Type_vector(
    /* count = */ A.numCols,
    /* blocklength = */ 1,
    /* stride = */ A.incCol,
    /* element type = */ get_type(A(0,0)),
    /* newly created type = */ &rowtype);

    /* test for rowmajor
       assuming that the matrix (view) is consecutively in memory
       that means no column gaps*/
    if (A.incRow == A.numCols) {
        MPI_Type_commit(&rowtype);
        return rowtype;
    } else {
        MPI_Datatype resize_rowtype;
        MPI_Type_create_resized(rowtype, 0,
                               A.incRow*sizeof(ElementType),
                               &resize_rowtype);
        MPI_Type_commit(&resize_rowtype);
        return resize_rowtype;
    }

}

template<typename Matrix>
typename std::enable_if<hpc::matvec::IsGeMatrix<Matrix>::value,
         MPI_Datatype>::type
get_type(const Matrix& A)
{
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
typename std::enable_if<hpc::matvec::IsGeMatrix<MA>::value &&
                        hpc::matvec::IsGeMatrix<MB>::value &&
                        std::is_same<typename MA::ElementType,
                        typename MB::ElementType>::value, int>::type
scatter_by_row(const MA& A, MB& B, int root, MPI_Comm comm)
{
    assert(A.numCols == B.numCols);

    int nof_processes;
    MPI_Comm_size(comm, &nof_processes);

    int rank;
    MPI_Comm_rank(comm, &rank);

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
    return MPI_Scatterv((void*) &A(0, 0), &counts[0], &offsets[0],
                        rowtype_A,
                        &B(0, 0), recvcount, rowtype_B,
                        root, comm);
}

template<typename MA, typename MB>
typename std::enable_if<hpc::matvec::IsGeMatrix<MA>::value &&
                        hpc::matvec::IsGeMatrix<MB>::value &&
                        std::is_same<typename MA::ElementType,
                        typename MB::ElementType>::value, int>::type
gather_by_row(const MA& A, MB& B, int root, MPI_Comm comm)
{
    assert(A.numCols == B.numCols);

    int nof_processes;
    MPI_Comm_size(comm, &nof_processes);

    int rank;
    MPI_Comm_rank(comm, &rank);

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
                       &B(0, 0), &counts[0], &offsets[0], rowtype_B,
                       root, comm);
}

} } // namespaces mpi, hpc

#endif // HPC_MPI_MATRIX_H
