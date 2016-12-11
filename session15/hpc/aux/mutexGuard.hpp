#ifndef INC_MUTEX_GUARD_HPP
#define INC_MUTEX_GUARD_HPP 1

#include <mutex>

class MutexGuard
{
    std::mutex& mutex;

    public:
        MutexGuard(std::mutex& mutex): mutex(mutex) {
            mutex.lock();
        }

        ~MutexGuard() {
            mutex.unlock();
        }

};

#endif      //INC_MUTEX_GUARD_HPP
