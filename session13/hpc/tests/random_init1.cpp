#include <random>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/print.h>
#include <hpc/matvec/apply.h>

template<typename MA>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
randomInit(MA& A) {
   using ElementType = typename MA::ElementType;
   using Index = typename MA::Index;

   std::random_device random;
   std::mt19937 mt(random());
   std::uniform_real_distribution<ElementType> uniform(-100,100);

   auto randFunc = [=] (ElementType &val,Index i, Index j) mutable
       -> void {
       val = uniform(mt);
   };
   hpc::matvec::apply(A, randFunc);
}

int main() {
   using namespace hpc::matvec;
   GeMatrix<double> A(7, 8);
   randomInit(A);
   print(A, "A");
}
