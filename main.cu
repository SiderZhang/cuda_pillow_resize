#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>

#include "libcuda.h"
#include "download/download.h"

#include <sys/timeb.h>
#include <unistd.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr<<"arguments: input_dir output_dir thread_count"<<std::endl;
        return -1;
    }
    char* input_dir = argv[1];
    char* output_dir = argv[2];
    std::string temp_path(output_dir);
    int thread_count = atoi(argv[3]);

    init(thread_count, temp_path);

    while (true) {
        timeb t;
        ftime(&t);
        std::vector<std::string> fileNames;
        read_dir_download(input_dir, fileNames);
        submit_download_job(fileNames);
        wait_downloading_over(fileNames);

        long t1 = t.time * 1000 + t.millitm;
        int count = 0;
        while(true) {
            ImageBase* image_base = get_next_image_base();
            if (image_base == nullptr)
                break;

//            std::string output_file_abs_path = std::string(output_dir) + "/" + image_base->image_suffix;

            int ret = image_process(image_base, image_base->image_path, 224, 224, 256, 256);
            if (ret != 0) {
                delete image_base;
                continue;
            }

            ret = remove(image_base->image_path.c_str());
            if (ret != 0) {
                std::cerr<<"failed to delete processed image file "<<image_base->image_path.c_str()<<std::endl;
                delete image_base;
                continue;
            }

            delete image_base;
            count++;
        }
        ftime(&t);
        long t2 = t.time * 1000 + t.millitm;

        if (count != 0) {
            std::cout << "process images " << count << " for time " << t2 - t1 << " millis at "<< t2 << std::endl;
        }
        usleep(10000);
    }

//    unsigned char* data;
//    read2("/home/siderzhang/file/9.jpg", &data);

//    if (argc < 1) {
//        std::cerr<<"arguments: input_filename"<<std::endl;
//        return -1;
//    }
//    const char* input_file = argv[1];
//    std::string output_filename_prefix = std::string("hello");
//    image_process(input_file, output_filename_prefix, 224, 224, 256, 256);
//    return 0;
}
