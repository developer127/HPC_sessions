#ifndef HPC_MPI_VECTOR_H
#define HPC_MPI_VECTOR_H 1

#include <mpi.h>
#include <hpc/matvec/densevector.h>
#include <hpc/matvec/isdensevector.h>
#include <hpc/mpi/fundamental.h>

namespace hpc { namespace mpi {

template<typename Vector>
typename std::enable_if<hpc::matvec::IsDenseVector<Vector>::value,
   MPI_Datatype>::type
get_type(const Vector& vector) {
   using ElementType = typename Vector::ElementType;
   MPI_Datatype datatype;
   MPI_Type_vector(
      /* count = */ vector.length,
      /* blocklength = */ 1,
      /* stride = */ vector.inc,
      /* element type = */ get_type(vector(0)),
      /* newly created type = */ &datatype);
   MPI_Type_commit(&datatype);
   return datatype;
}

} } // namespaces mpi, hpc

#endif // HPC_MPI_VECTOR_H
