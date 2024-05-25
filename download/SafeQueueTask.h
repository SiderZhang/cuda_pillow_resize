//
// Created by siderzhangPC on 2024/5/23.
//

#ifndef CURLDOWNLOAD_SAFEQUEUETASK_H
#define CURLDOWNLOAD_SAFEQUEUETASK_H

#include <queue>
#include <mutex>

#include <mutex>
#include <queue>
#include <functional>
#include <future>
#include <thread>
#include <utility>
#include <vector>

#include "Task.h"

// Thread safe implementation of a Queue using a std::queue
class SafeQueueTask
{
private:
    std::queue<Task> m_queue; //利用模板函数构造队列
    
    std::mutex m_mutex; // 访问互斥信号量
    
public:
    SafeQueueTask() {}
    SafeQueueTask(SafeQueueTask &&other) {}
    ~SafeQueueTask() {}
    
    bool empty() // 返回队列是否为空
    {
        std::unique_lock<std::mutex> lock(m_mutex); // 互斥信号变量加锁，防止m_queue被改变
        
        return m_queue.empty();
    }
    
    int size()
    {
        std::unique_lock<std::mutex> lock(m_mutex); // 互斥信号变量加锁，防止m_queue被改变
        
        return m_queue.size();
    }
    
    // 队列添加元素
    void enqueue(Task& t)
    {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_queue.push(t);
    }
    
    // 队列取出元素
    Task dequeue()
    {
        std::unique_lock<std::mutex> lock(m_mutex); // 队列加锁
        Task t = m_queue.front(); // 取出队首元素，返回队首元素值，并进行右值引用
        m_queue.pop(); // 弹出入队的第一个元素
        return t;
    }
};

#endif //CURLDOWNLOAD_SAFEQUEUE_H
