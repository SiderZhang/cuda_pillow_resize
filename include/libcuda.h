//
// Created by siderzhangPC on 2024/4/22.
//

#ifndef LIBJPEGTEST_LIBCUDA_H
#define LIBJPEGTEST_LIBCUDA_H

float* resize(const char* filename, unsigned int target_xsize, unsigned int target_ysize);
void read2(const char* filepath, unsigned char **output_buffer, unsigned int& width, unsigned int& height, unsigned int& channels);
void onImageRead(unsigned char* input, unsigned char **output_buffer, unsigned int width, unsigned int height, unsigned int channels);

#endif //LIBJPEGTEST_LIBCUDA_H
