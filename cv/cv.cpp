//
// Created by siderzhangPC on 2024/5/12.
//

#include "opencv4/opencv2/opencv.hpp"
#include <string>
#include "libcuda.h"
#include <iostream>

using namespace cv;
using namespace std;

void read2(const char* filepath, unsigned char **output_buffer, unsigned int& width, unsigned int& height, unsigned int& channels){
    std::string path(filepath);
//    cv::Mat input = cv::imread(path, cv::IMREAD_COLOR);
    cv::Mat src_input = cv::imread(path);
    cv::Mat input;
    cvtColor(src_input, input, cv::COLOR_BGR2RGB);
//    cv::Mat input = cv2.cvtColor(src_input, cv2.COLOR_BGR2RGB).clone();
    cout<<"channels: "<<input.channels()<<"; type: "<<input.type()<<"; size: "<<input.size()<<"; step: "<<input.step[0]<<endl;
    width = input.cols;
    height = input.rows;
    channels = input.channels();
    onImageRead(input.ptr(), output_buffer, width, height, channels);
    input.release();
//    src_input.release();
    printf("Hello\n");
}