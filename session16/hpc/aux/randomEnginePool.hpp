#ifndef INC_RANDOM_ENGINE_POOL
#define INC_RANDOM_ENGINE_POOL 1

#include <list>

namespace hpc { namespace aux {

template<typename T>
struct RandomEnginePool {
    using EngineType = T;
    T get() {
        /* check if we have a free engine in unused */
        {
        std::lock_guard<std::mutex> lock(mutex);
        if (unused.size() > 0) {
            T rg = std::move(unused.front());
            unused.pop_front();
            return rg;
        }
        }
        /* prepare new random generator */
        return T(r());
    }
    void free(T&& engine) {                       //asumes rvalue reference
        std::lock_guard<std::mutex> lock(mutex);
        unused.push_back(engine);
    }
    private:
        std::mutex mutex;
        std::random_device r;
        std::list<T> unused;
};
}}      // namespace hpc::aux
#endif  // INC_RANDOM_ENGINE_POOL
