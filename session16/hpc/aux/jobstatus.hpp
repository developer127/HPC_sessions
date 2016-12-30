#ifndef INC_AUX_JOBSTATUS_HPP
#define INC_AUX_JOBSTATUS_HPP 1

namespace hpc { namespace aux {

struct JobStatus {
    public:
        JobStatus() : finished(false) {
        }
        void is_finished() {
            std::unique_lock<std::mutex> lock(mutex);
            assert(!finished);
            finished = true;
            cv.notify_all();
        }
        void wait() {
            std::unique_lock<std::mutex> lock(mutex);
            if (!finished) {
                cv.wait(lock);
            }
        }
    private:
        std::mutex mutex;
        std::condition_variable cv;
        bool finished;
};

} }         // namespace hpc::aux

#endif      //INC_AUX_JOBSTATUS_HPP
