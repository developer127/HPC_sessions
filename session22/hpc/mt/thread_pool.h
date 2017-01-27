/* 
   Copyright (c) 2015, 2016 Andreas F. Borchert
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
   FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
   COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
   BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

/*
   This header-only C++ package provides hpc::mt::thread_pool which
   assigns incoming tasks to a fixed number of threads.

   Following methods are supported:

    * ThreadPool(unsigned int nofthreads)
      construct thread pool with the given number of threads

    * ~ThreadPool()
      the destructor waits for all threads to finish

    * unsigned int get_num_threads()
      return the number of threads

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
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace hpc { namespace mt {

   class ThreadPool {
      public:
	 ThreadPool(unsigned int nofthreads) :
	       threads(nofthreads), active(0),
	       joining(false), terminating(false) {
	    for (auto& t: threads) {
	       t = std::thread([=]() mutable -> void {
		  for (;;) {
		     std::function<void()> task;
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

	 unsigned int get_num_threads() const {
	    return threads.size();
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
	 auto submit(F&& task_function)
	       -> std::future<decltype(task_function())> {
	    using T = decltype(task_function());
	    /* a std::function object cannot be constructed from
	       a std::packaged_task object as std::function requires
	       the callable object to be copy-constructible;
	       packaged tasks are, however, just move-constructible;
	       the following workaround puts the packaged task
	       behind a shared pointer and passes a simple lambda
	       object to the std::function constructor
	    */
	    auto task = std::make_shared<std::packaged_task<T()>>(
	       std::forward<F>(task_function));
	    std::future<T> result = task->get_future();
	    std::lock_guard<std::mutex> lock(mutex);
	    tasks.push_back([task]() { (*task)(); } );
	    cv.notify_one();
	    return result;
	 }

      private:
	 std::mutex mutex;
	 std::condition_variable cv;
	 std::vector<std::thread> threads; // fixed number of threads
	 std::list<std::function<void()>> tasks; // submitted tasks
	 unsigned int active; // number of threads that executing a task
	 bool joining; // initially false, set to true if join() is invoked
	 bool terminating; // initially false, set to true by terminate()
   };

} } // namespaces mt and hpc

#endif
