#ifndef HPC_MPI_RECEIVE_VEC_HPP
#define HPC_MPI_RECEICE_VEC_HPP 1

#include <hpc/mpi/get_type.hpp>
#include <cassert>

namespace hpc { namespace mpi {

template<typename Vector>
void
receiveVec(Vector& vector, int source, int tag)
{
    assert(source >= 0);
    MPI_Status status;
    MPI_Datatype datatype = get_type(vector);
    MPI_Recv(&vector(0), 1, datatype,
             source, tag, MPI_COMM_WORLD, &status);
}

} }         // namespace hpc::mpi

#endif      //HPC_MPI_RECEIVE_VEC_HPP
