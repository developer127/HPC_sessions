#ifndef HPC_MPI_VECTOR_HPP
#define HPC_MPI_VECTOR_HPP 1

#include <cassert>
#include <hpc/mpi/fundamental.h>
#include <hpc/matvec/isdensevector.h>

namespace hpc { namespace mpi {

template<typename Vector>
typename std::enable_if<hpc::matvec::IsDenseVector<Vector>::value,
    MPI_Datatype>::type
get_type(const Vector& vector)
{
    using ElementType = typename Vector::ElementType;
    MPI_Datatype datatype;
    MPI_Type_vector(vector.length, 1, vector.inc,
                    fundamental_type<ElementType>().get(),
                    &datatype);
    MPI_Type_commit(&datatype);
    return datatype;
}


template<typename Vector>
void
sendVec(const Vector& vector, int dest, int tag)
{
    assert(dest >= 0);
    MPI_Datatype datatype = get_type(vector);
    MPI_Send(&vector(0), 1, datatype,
             dest, tag, MPI_COMM_WORLD);
}


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

#endif      //HPC_MPI_VECTOR_HPP
