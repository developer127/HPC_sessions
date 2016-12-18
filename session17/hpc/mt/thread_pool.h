/*
   This header-only C++ package provides hpc::mt::thread_pool which
   assigns incoming tasks to a fixed number of threads.  The template
   parameter determines the return type of a task.

   Following methods are supported:

    * ThreadPool(unsigned int nofthreads)
      construct thread pool with the given number of threads

    * ~ThreadPool()
      the destructor waits for all threads to finish

    * void join()
      wait for all threads to finish

    * void terminate()
      request a speedy termination, i.e. submitted but not yet assigned
      tasks remain undone; threads that wait for the corresponding futures
      will see broken promises

    * std::future<T> submit(F task_function)
      submit a task which is to be executed by one of the threads
      of the pool; the future objects allows to synchronize with
      the completion of the task and to receive the return value
      of the submitted task
*/

#ifndef HPC_MT_THREAD_POOL_H
#define HPC_MT_THREAD_POOL_H 1

#include <condition_variable>
#include <future>
#include <list>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace hpc { namespace mt {

   template<typename T>
   class ThreadPool {
      public:
	 ThreadPool(unsigned int nofthreads) :
	       threads(nofthreads), active(0),
	       joining(false), terminating(false) {
	    for (auto& t: threads) {
	       t = std::thread([=]() mutable -> void {
		  for (;;) {
		     std::packaged_task<T()> task;
		     /* fetch next task, if there is any */
		     {
			std::unique_lock<std::mutex> lock(mutex);
			while (!terminating &&
			      tasks.empty() &&
			      (active > 0 || !joining)) {
			   cv.wait(lock);
			}
			if (terminating || tasks.empty()) {
			   terminating = true;
			   break; /* all done */
			}
			task = std::move(tasks.front());
			tasks.pop_front();
			++active;
		     }
		     /* process task */
		     task();
		     /* decrement number of active threads */
		     {
			std::unique_lock<std::mutex> lock(mutex);
			--active;
		     }
		  }
		  cv.notify_all();
	       });
	    }
	 }
	 ~ThreadPool() {
	    join();
	 }

	 void join() {
	    if (!joining && !terminating) {
	       std::unique_lock<std::mutex> lock(mutex);
	       joining = true;
	    }
	    cv.notify_all();
	    for (auto& t: threads) if (t.joinable()) t.join();
	 }

	 void terminate() {
	    {
	       std::unique_lock<std::mutex> lock(mutex);
	       terminating = true;
	       /* we make sure that all promises left are considered broken
		  by emptying the list of remaining tasks;
		  if we do not do it now, the waiting threads would
		  have to wait until this object is destructed
	       */
	       tasks.empty();
	    }
	    cv.notify_all();
	 }

	 template<typename F>
	 std::future<T> submit(F task_function) {
	    std::packaged_task<T()> task(task_function);
	    std::future<T> result = task.get_future();
	    std::lock_guard<std::mutex> lock(mutex);
	    tasks.push_back(std::move(task));
	    cv.notify_one();
	    return result;
	 }
      private:
	 std::vector<std::thread> threads;
	 std::mutex mutex;
	 std::list<std::packaged_task<T()>> tasks;
	 std::condition_variable cv;
	 unsigned int active;
	 bool joining;
	 bool terminating;
   };

} } // namespaces mt and hpc

#endif
