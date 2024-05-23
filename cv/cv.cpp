//
// Created by siderzhangPC on 2024/5/12.
//

#include "opencv2/opencv.hpp"
#include <string>
#include "libcuda.h"
#include <iostream>

using namespace cv;
using namespace std;

int load_image_file(const char* filepath, unsigned char **output_buffer, unsigned int& width, unsigned int& height, unsigned int& channels){
    std::string path(filepath);
//    cv::Mat input = cv::imread(path, cv::IMREAD_COLOR);
    cv::Mat src_input = cv::imread(path);
    if (src_input.empty()) {
        std::cout<<"Empty content got! Invalid image file: "<< path << std::endl;
        return -1;
    }

    cv::Mat input;
    cvtColor(src_input, input, cv::COLOR_BGR2RGB);
    width = input.cols;
    height = input.rows;
    channels = input.channels();
    onImageRead(input.ptr(), output_buffer, width, height, channels);
    input.release();

    return 0;
}