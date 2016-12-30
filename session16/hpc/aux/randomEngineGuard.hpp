#ifndef INC_RANDOM_ENGINE_GUARD
#define INC_RANDOM_ENGINE_GUARD 1

#include <hpc/aux/randomEnginePool.hpp>

namespace hpc { namespace aux {

template<typename T>
struct RandomEngineGuard {
    using EngineType = T;

    RandomEngineGuard(RandomEnginePool<EngineType>& pool):
            pool(pool),
            engine(pool.get()) {}

    EngineType& get() {
        return engine;
    }

    ~RandomEngineGuard()
    {
        pool.free(std::move(engine));
    }
    private:
        RandomEnginePool<EngineType>& pool;
        EngineType engine;
};
}}      // namespace hpc::aux
#endif  // INC_RANDOM_ENGINE_GUARD
