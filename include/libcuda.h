//
// Created by siderzhangPC on 2024/4/22.
//

#ifndef LIBJPEGTEST_LIBCUDA_H
#define LIBJPEGTEST_LIBCUDA_H

extern "C" {
    float* resize(const char* filename, unsigned int target_xsize, unsigned int target_ysize);

};

#endif //LIBJPEGTEST_LIBCUDA_H
