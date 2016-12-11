#include <random>
#include <hpc/matvec/gematrix.h>
#include <hpc/matvec/apply.h>

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

