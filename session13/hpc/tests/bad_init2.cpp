#include <thread>
#include <random>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/apply.h>
#include <hpc/matvec/print.h>

template<typename MA, typename Func>
typename std::enable_if<hpc::matvec::IsRealGeMatrix<MA>::value, void>::type
randomInit(MA& A, Func& f) {
   using ElementType = typename MA::ElementType;
   using Index = typename MA::Index;
   
   hpc::matvec::apply(A, [&](ElementType& val, Index i, Index j) -> void {
      val = f();
   });
}

int main() {
   std::random_device random;
   std::mt19937 mt(random());
   std::uniform_real_distribution<double> uniform(-100,100);
   auto rand = [=]() mutable -> double { return uniform(mt); };

   using namespace hpc::matvec;
   GeMatrix<double> A(1000, 1000);
   GeMatrix<double> B(1000, 1000);
   GeMatrix<double> C(1000, 1000);

   /* start three threads that initialize A, B, and C */
   std::thread t1([&](){ randomInit(A, rand); });
   std::thread t2([&](){ randomInit(B, rand); });
   std::thread t3([&](){ randomInit(C, rand); });

   /* wait until they are finished */
   t1.join(); t2.join(); t3.join();

   /* print a small block of each of the initialized matrices */
   auto A1 = A(500, 500, 5, 5);
   auto B1 = B(500, 500, 5, 5);
   auto C1 = C(500, 500, 5, 5);
   print(A1, "A1");
   print(B1, "B1");
   print(C1, "C1");
}
