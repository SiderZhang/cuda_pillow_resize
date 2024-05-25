//
// Created by siderzhangPC on 2024/4/22.
//

#ifndef LIBJPEGTEST_LIBCUDA_H
#define LIBJPEGTEST_LIBCUDA_H

//#include

struct ImageBase {
    std::string image_path;
    std::string image_suffix;
    unsigned char* pixel_data_d = nullptr;
    unsigned int source_xsize = 0;
    unsigned int source_ysize = 0;
    unsigned int channels = 0;
};

float* resize(const char* filename, unsigned int target_xsize, unsigned int target_ysize);
//int load_image_file(const char* filepath, unsigned char **output_buffer, unsigned int& width, unsigned int& height, unsigned int& channels);
int load_image_file(const char* filepath, ImageBase* image_base);
void onImageRead(unsigned char* input, unsigned char **output_buffer, unsigned int width, unsigned int height, unsigned int channels);
void readDir(const char* dirPath, std::vector<std::string>& filenames, std::string output_path);
void read_dir_download(std::string dirPath, std::vector<std::string>& filenames);
int image_process(ImageBase* image_base, std::string& output_filename_prefix, unsigned int vit_size, unsigned int vis_size, unsigned int deit_size, unsigned int corp_size);
//int image_process(unsigned char* input, unsigned int source_xsize, unsigned int source_ysize, unsigned int channels, std::string& output_filename_prefix, unsigned int vit_size, unsigned int vis_size, unsigned int deit_size, unsigned int corp_size);
void onInvalidImage(unsigned char **output_buffer, unsigned int& width, unsigned int& height, unsigned int& channels);
void to_npy(std::string& outputFilename, float* array, unsigned int size, unsigned int channels);

#endif //LIBJPEGTEST_LIBCUDA_H
