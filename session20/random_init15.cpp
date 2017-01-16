#include <hpc/aux/repool.h>
#include <hpc/aux/slices.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/asum.h>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/print.h>
#include <omp.h>

template<typename MA, typename Pool>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
randomInit(MA& A, Pool& pool) {
   using ElementType = typename MA::ElementType;
   using Index = typename MA::Index;
   using EngineType = typename Pool::EngineType;

   std::uniform_real_distribution<double> uniform(-100, 100);
   hpc::aux::RandomEngineGuard<EngineType> guard(pool);
   auto& engine(guard.get());

   hpc::matvec::apply(A, [&](ElementType& val, Index i, Index j) -> void {
      val = uniform(engine);
   });
}

template<typename MA, typename RandomEnginePool>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
mt_randomInit(MA& A, RandomEnginePool& repool) {
   using ElementType = typename MA::ElementType;
   using Index = typename MA::Index;
   using namespace hpc::aux;

   #pragma omp parallel for
   for (Index row = 0; row < A.numRows; ++row) {
      auto A_ = A(row, 0, 1, A.numCols);
      randomInit(A_, repool);
   }
}

int main() {
   using namespace hpc::matvec;
   using namespace hpc::aux;
   RandomEnginePool<std::mt19937> repool(4);

   GeMatrix<double> A(51, 7);
   mt_randomInit(A, repool);
}
