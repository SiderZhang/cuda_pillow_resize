//
// Created by siderzhangPC on 2024/4/22.
//

#ifndef LIBJPEGTEST_LIBCUDA_H
#define LIBJPEGTEST_LIBCUDA_H


float* resize(const char* filename, unsigned int target_xsize, unsigned int target_ysize);
int load_image_file(const char* filepath, unsigned char **output_buffer, unsigned int& width, unsigned int& height, unsigned int& channels);
void onImageRead(unsigned char* input, unsigned char **output_buffer, unsigned int width, unsigned int height, unsigned int channels);
void readDir(const char* dirPath, std::vector<std::string>& filenames, std::string output_path);

#endif //LIBJPEGTEST_LIBCUDA_H
