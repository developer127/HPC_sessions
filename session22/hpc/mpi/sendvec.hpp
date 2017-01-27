#ifndef HPC_MPI_SEND_VEC_HPP
#define HPC_MPI_SEND_VEC_HPP 1

#include <hpc/mpi/get_type.hpp>
#include <cassert>

namespace hpc { namespace mpi {

template<typename Vector>
void
sendVec(const Vector& vector, int dest, int tag)
{
    assert(dest >= 0);
    MPI_Datatype datatype = get_type(vector);
    MPI_Send(&vector(0), 1, datatype,
             dest, tag, MPI_COMM_WORLD);
}

} }         // namespace hpc::mpi

#endif      //HPC_MPI_SEND_VEC_HPP
