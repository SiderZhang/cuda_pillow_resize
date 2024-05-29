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
typedef void (*Func)(std::string&, std::string output_dir, SafeQueue*);

class Task {
private :
    std::string image_suffix;
    std::string output_dir;
    Func f;
    SafeQueue* queue;

public:
    Task(): f(nullptr) {

    }

    Task(std::string _image_suffix, std::string _output_dir, Func _f, SafeQueue *_queue) : image_suffix(std::move(_image_suffix)), output_dir(_output_dir), f(_f), queue(_queue) {

    }


    void operator()() {
        f(image_suffix, output_dir, queue);
    }

    std::string getData() {
        return image_suffix;
    }
};

#endif //CURLDOWNLOAD_TASK_H
