#include <cassert>
#include <condition_variable>
#include <functional>
#include <list>
#include <mutex>
#include <thread>
#include <utility>
#include <fmt/printf.hpp>

struct Worker {
    public:
        using Job = std::function<void()>;
        Worker(): isfinished(false) {
            t = std::thread([=]() { process_jobs(); });
        }
        ~Worker() {
            isfinished = true;
            t.join();
        }
        void submit(Job job) {
            std::unique_lock<std::mutex> lock(mutex);
            jobs.push_back(std::move(job));
            cv.notify_one();
        }
    private:
        std::thread t;
        std::mutex mutex;
        std::condition_variable cv;
        std::list<Job> jobs;
        bool isfinished;

        void process_jobs() {
            for(;;) {
                Job job;
                /* fetch job */
                {
                    std::unique_lock<std::mutex> lock(mutex);
                    if (jobs.empty() && isfinished == false) {
                        cv.wait(lock);
                    }
                    if (jobs.empty() && isfinished == true) {
                        break;
                    }
                    job = std::move(jobs.front());
                    jobs.pop_front();
                }
                /* execute job */
                job();
            }
        }
};

int main() {
    Worker worker;
    worker.submit([]() { fmt::printf("Hi, this is your first job!\n"); });
    worker.submit([]() { fmt::printf("Now you got another job.\n"); });
}
