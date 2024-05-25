//
// Created by siderzhangPC on 2024/5/12.
//

#include "opencv2/opencv.hpp"
#include <string>
#include "libcuda.h"
#include <iostream>

using namespace cv;
using namespace std;

int load_image_file(const char* filepath, ImageBase* image_base) {// unsigned char **output_buffer, unsigned int& width, unsigned int& height, unsigned int& channels){
    std::string path(filepath);
//    cv::Mat input = cv::imread(path, cv::IMREAD_COLOR);
    cv::Mat src_input = cv::imread(path);
    if (src_input.empty()) {
        std::cout<<"Empty content got! Invalid image file: "<< path << std::endl;
        return -1;
    }

    cv::Mat input;
    cvtColor(src_input, input, cv::COLOR_BGR2RGB);
    image_base->source_xsize = input.cols;
    image_base->source_ysize = input.rows;
    image_base->channels = input.channels();
    onImageRead(input.ptr(), &(image_base->pixel_data_d), image_base->source_xsize, image_base->source_ysize, image_base->channels);
    input.release();

    return 0;
}