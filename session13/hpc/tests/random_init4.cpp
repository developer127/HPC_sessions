#include <thread>
#include <random>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/print.h>

template<typename MA>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
randomInit(MA& A) {
   using ElementType = typename MA::ElementType;
   using Index = typename MA::Index;

   std::random_device random;
   std::mt19937 mt(random());
   std::uniform_real_distribution<ElementType> uniform(-100,100);
   
   hpc::matvec::apply(A, [&](ElementType& val, Index i, Index j) -> void {
      val = uniform(mt);
   });
}

int main() {
   using namespace hpc::matvec;
   GeMatrix<double> A(1000, 1000);

   randomInit(A);

   /* print a small block of each of the initialized matrices */
   auto A1 = A(500, 500, 5, 5);
   print(A1, "A1");
}
