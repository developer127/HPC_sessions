#ifndef HPC_MPI_GET_TYPE_HPP
#define HPC_MPI_GET_TYPE_HPP 1

#include <hpc/mpi/fundamental.h>

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

} }         // namespace hpc::mpi

#endif      //HPC_MPI_GET_TYPE_HPP
