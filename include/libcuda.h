//
// Created by siderzhangPC on 2024/4/22.
//

#ifndef LIBJPEGTEST_LIBCUDA_H
#define LIBJPEGTEST_LIBCUDA_H

extern "C" {
    void resize(unsigned int source_xsize, unsigned int source_ysize, unsigned char* image_pixel_data, int channels,
            unsigned int target_xsize, unsigned int target_ysize, float *normalized_pixel_data);

};

#endif //LIBJPEGTEST_LIBCUDA_H
