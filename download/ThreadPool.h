//
// Created by siderzhangPC on 2024/5/23.
//

#ifndef CURLDOWNLOAD_THREADPOOL_H
#define CURLDOWNLOAD_THREADPOOL_H

#include "SafeQueue.h"
#include "SafeQueueTask.h"
#include "Task.h"

class ThreadPool
{
private:
    class ThreadWorker // 内置线程工作类
    {
    private:
        int m_id; // 工作id

        ThreadPool *m_pool; // 所属线程池
    public:
        // 构造函数
        ThreadWorker(ThreadPool *pool, const int id) : m_pool(pool), m_id(id)
        {
        }

        // 重载()操作
        void operator()()
        {
            Task func; // 定义基础函数类func

            bool dequeued; // 是否正在取出队列中元素

            while (!m_pool->m_shutdown)
            {
                {
                    // 为线程环境加锁，互访问工作线程的休眠和唤醒
                    std::unique_lock<std::mutex> lock(m_pool->m_conditional_mutex);

                    // 如果任务队列为空，阻塞当前线程
                    if (m_pool->m_queue.empty())
                    {
                        m_pool->m_conditional_lock.wait(lock); // 等待条件变量通知，开启线程
                    }

                    // 取出任务队列中的元素
                    func = m_pool->m_queue.dequeue();
                }

                // 如果成功取出，执行工作函数
                func();
            }
        }
    };

    bool m_shutdown; // 线程池是否关闭

    SafeQueueTask m_queue; // 执行函数安全队列，即任务队列

    std::vector<std::thread> m_threads; // 工作线程队列

    std::mutex m_conditional_mutex; // 线程休眠锁互斥变量

    std::condition_variable m_conditional_lock; // 线程环境锁，可以让线程处于休眠或者唤醒状态

public:
    // 线程池构造函数
    ThreadPool(const int n_threads = 4)
            : m_threads(std::vector<std::thread>(n_threads)), m_shutdown(false)
    {
    }

    ThreadPool(const ThreadPool &) = delete;

    ThreadPool(ThreadPool &&) = delete;

    ThreadPool &operator=(const ThreadPool &) = delete;

    ThreadPool &operator=(ThreadPool &&) = delete;

    // Inits thread pool
    void init()
    {
        for (int i = 0; i < m_threads.size(); ++i)
        {
            m_threads.at(i) = std::thread(ThreadWorker(this, i)); // 分配工作线程
        }
    }

    // Waits until threads finish their current task and shutdowns the pool
    void shutdown()
    {
        m_shutdown = true;
        m_conditional_lock.notify_all(); // 通知，唤醒所有工作线程

        for (int i = 0; i < m_threads.size(); ++i)
        {
            if (m_threads.at(i).joinable()) // 判断线程是否在等待
            {
                m_threads.at(i).join(); // 将线程加入到等待队列
            }
        }
    }

    // Submit a function to be executed asynchronously by the pool
    void submit(Task& task)
    {
        // 队列通用安全封包函数，并压入安全队列
        m_queue.enqueue(task);

//        std::cout<<"submitted" << task->getData()<<std::endl;

        // 唤醒一个等待中的线程
        m_conditional_lock.notify_one();
    }
};

#endif //CURLDOWNLOAD_THREADPOOL_H
