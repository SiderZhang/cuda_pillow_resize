//
// Created by siderzhangPC on 2024/5/23.
//

#ifndef CURLDOWNLOAD_TASK_H
#define CURLDOWNLOAD_TASK_H


#include <string>
#include <utility>
#include "SafeQueue.h"
#include "libcuda.h"

//typedef ImageBase* (*Func)(std::string&);
typedef void (*Func)(std::string&, SafeQueue*);

class Task {
private :
    std::string image_suffix;
    Func f;
    SafeQueue* queue;
public:
    Task(): f(nullptr) {

    }

    Task(std::string _image_suffix, Func _f, SafeQueue *_queue) : image_suffix(std::move(_image_suffix)), f(_f), queue(_queue) {

    }


    void operator()() {
        f(image_suffix, queue);
    }

    std::string getData() {
        return image_suffix;
    }
};

#endif //CURLDOWNLOAD_TASK_H
