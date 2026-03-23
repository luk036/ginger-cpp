// ThreadPool.h
#pragma once

#include <condition_variable>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

/**
 * @brief ThreadPool
 *
 * A simple thread pool implementation that allows tasks to be enqueued and
 * executed asynchronously by a fixed number of worker threads.
 *
 * The pool creates a specified number of worker threads at construction.
 * Tasks can be enqueued using the enqueue() method.
 *
 * @note The pool will automatically shut down when destroyed, waiting for
 * all pending tasks to complete.
 */
class ThreadPool {
  public:
    /**
     * @brief Construct a new ThreadPool
     *
     * Creates a thread pool with the specified number of worker threads.
     *
     * @param threads The number of worker threads to create
     */
    ThreadPool(size_t threads);
    
    /**
     * @brief Destroy the ThreadPool
     *
     * Stops the pool and waits for all worker threads to complete their current tasks.
     */
    ~ThreadPool();

    /**
     * @brief Enqueue a new task to the thread pool
     *
     * Adds a new work item to the pool's task queue. The task will be executed
     * by one of the worker threads.
     *
     * @tparam F The callable type of the task
     * @tparam Args The argument types for the task
     * @param f The callable to execute
     * @param args The arguments to pass to the callable
     */
    template <class F, class... Args> void enqueue(F&& f, Args&&... args);

  private:
    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> tasks_;
    std::mutex queue_mutex_;
    std::condition_variable condition_;
    bool stop_;
};

// ThreadPool.cpp
inline ThreadPool::ThreadPool(size_t threads) : stop_(false) {
    for (size_t i = 0; i < threads; ++i) {
        workers_.emplace_back([this]() {
            while (true) {
                std::function<void()> task;
                {
                    std::unique_lock<std::mutex> lock(this->queue_mutex_);
                    this->condition_.wait(lock,
                                          [this] { return this->stop_ || !this->tasks_.empty(); });
                    if (this->stop_ && this->tasks_.empty()) return;
                    task = std::move(this->tasks_.front());
                    this->tasks_.pop();
                }
                task();
            }
        });
    }
}

inline ThreadPool::~ThreadPool() {
    {
        std::unique_lock<std::mutex> lock(queue_mutex_);
        stop_ = true;
    }
    condition_.notify_all();
    for (std::thread& worker : workers_) {
        worker.join();
    }
}

template <class F, class... Args> inline void ThreadPool::enqueue(F&& f, Args&&... args) {
    auto task = std::bind(std::forward<F>(f), std::forward<Args>(args)...);
    {
        std::lock_guard<std::mutex> lock(queue_mutex_);
        tasks_.emplace(task);
    }
    condition_.notify_one();
}
