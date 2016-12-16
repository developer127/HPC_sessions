#ifndef INC_RANDOM_ENGINE_GUARD_FIX_SIZE_HPP
#define INC_RANDOM_ENGINE_GUARD_FIX_SIZE_HPP 1

#include <hpc/aux/randomEnginePoolFixSize.hpp>

namespace hpc { namespace aux {

template<typename T>
struct RandomEngineGuard {
    using EngineType = T;

    RandomEngineGuard(RandomEnginePool<EngineType>& pool):
            pool(pool),
            engine(pool.get()) {}

    T& get() {
        return engine;
    }

    ~RandomEngineGuard()
    {
        pool.free(engine);
    }
    private:
        RandomEnginePool<T>& pool;
        T& engine;
};
}}      // namespace hpc::aux
#endif  // INC_RANDOM_ENGINE_GUARD_FIX_SIZE_HPP
