#ifndef INC_AUX_TPOOL_HPP
#define INC_AUX_TPOOL_HPP 1

#include <cassert>
#include <condition_variable>
#include <functional>
#include <list>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>
#include <fmt/printf.hpp>

namespace hpc {  namespace aux {

struct ThreadPool {
    public:
        using Job = std::function<void()>;

        ThreadPool(unsigned int nof_threads) :
            nof_threads(nof_threads), nof_active_threads(0),
            finished(false), threads(nof_threads) {
            for (auto& t: threads) {
                t = std::thread([=]() { process_jobs(); });
            }
        }

        ~ThreadPool() {
            {
                std::unique_lock<std::mutex> lock(mutex);
                finished = true;
            }
            cv.notify_all();
            for (auto& t: threads) {
                t.join();
            }
        }

        void submit(Job job) {
            std::unique_lock<std::mutex> lock(mutex);
            jobs.push_back(std::move(job));
            cv.notify_one();
        }
        unsigned int get_nof_threads() {
            std::unique_lock<std::mutex> lock(mutex);
            return nof_threads;
        }
    private:
        unsigned int nof_threads;
        unsigned int nof_active_threads;
        bool finished;
        std::vector<std::thread> threads;
        std::mutex mutex;
        std::condition_variable cv;
        std::list<Job> jobs;

        void process_jobs() {
            for(;;) {
                Job job;
                /* fetch job */
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    while (jobs.empty() && (!finished||nof_active_threads>0)){
                        cv.wait(lock);
                    }
                    if (jobs.empty() && nof_active_threads == 0 && finished) {
                        break;
                    }
                    job = std::move(jobs.front());
                    jobs.pop_front();
                    ++ nof_active_threads;
                }
                /* execute job */
                job();
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    -- nof_active_threads;
                }
            }
            /*if one thread finishes all finish*/
            cv.notify_all();
        }
};

} }         // namespace hpc::aux
#endif      //INC_AUX_TPOOL_HPP
