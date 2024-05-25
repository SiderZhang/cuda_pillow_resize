#include <string>
#include <cstring>
#include <vector>
#include <iostream>
#include <unistd.h>

#include <algorithm>
#include <dirent.h>
#include "sys/stat.h"
#include "npy.h"

using namespace std;

void scanDir(const string& base_path, const string& dir_path, std::vector<std::string>& filenames, std::string output_path) {
    DIR *pDir;
    struct dirent* ptr;
    std::string src_folder = std::string(base_path).append(dir_path);
    for (int i = 0;i < 5;i++) {
        pDir = opendir(src_folder.c_str());
        if (pDir != NULL) {
            break;
        }

        std::cout<<"Failed to open folder <" << src_folder << ">! "<<"Retry "<< (i + 1) << std::endl;
        usleep(5000);
    }

    std::string path_to_create = std::string(output_path).append(dir_path);
    if (access(path_to_create.c_str(), F_OK) != 0) {
        int ret = mkdir(path_to_create.c_str(), 0755);
        if (ret != 0) {
            if (ret != EEXIST) {
                std::cerr<<"failed to make directory for file: " << std::string(output_path) << ";" << std::string(dir_path) << ";" << ret <<std::endl;
            }
        }
    }

    struct stat s_buff;

    while((ptr = readdir(pDir)) != nullptr) {
        std::string fileName(ptr->d_name);
        std::string abs_path = std::string(base_path).append("/").append(dir_path).append("/").append(fileName);
        std::string relative_path = std::string(dir_path).append("/").append(fileName);

        std::string extension = fileName.substr(fileName.find_last_of(".") + 1);
        transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

        // ignore temp files
        if (extension.compare("tmp") == 0) {
            continue;
        }

        if (strcmp(ptr->d_name, ".") == 0 || strcmp(ptr->d_name, "..") == 0){
            continue;
        }

        stat(abs_path.c_str(), &s_buff);
        if (S_ISREG(s_buff.st_mode)) {
            filenames.push_back(relative_path);

            continue;
        }

        if (S_ISDIR(s_buff.st_mode)) {
            scanDir(base_path, relative_path, filenames, output_path);
            continue;
        }
    }
    closedir(pDir);
}

void readDir(const char* dirPath, std::vector<std::string>& filenames, std::string output_path) {
    scanDir(std::string(dirPath), std::string(""), filenames, output_path);

//    DIR *pDir;
//    struct dirent* ptr;
//    if(!(pDir = opendir(dirPath))){
//        std::cout<<"Folder doesn't Exist!"<<std::endl;
//        return;
//    }
//    struct stat s_buff;
//
//    while((ptr = readdir(pDir)) != nullptr) {
//        std::string fileName(ptr->d_name);
//        std::string path = std::string(dirPath) + "/" + fileName;
//
//        std::string extension = fileName.substr(fileName.find_last_of(".") + 1);
//        transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
//        printf("%s", ptr->d_name);
//        // ignore temp files
//        if (!extension.compare("tmp") == 0) {
//            continue;
//        }
//
//        if (strcmp(ptr->d_name, ".") == 0 || strcmp(ptr->d_name, "..") == 0){
//            continue;
//        }
//
//        stat(path.c_str(), &s_buff);
//        if (!S_ISREG(s_buff.st_mode)) {
//            continue;
//        }
//
//        filenames.push_back(fileName);
//    }
//    closedir(pDir);
}

void read_dir_download(std::string dirPath, std::vector<std::string>& filenames) {

    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(dirPath.c_str()))){
        std::cout<<"Folder doesn't Exist!"<<std::endl;
        return;
    }
    struct stat s_buff;

    while((ptr = readdir(pDir)) != nullptr) {
        std::string fileName(ptr->d_name);
        std::string tmp_dir_path = dirPath;
        std::string path = tmp_dir_path.append("/").append(fileName);

        std::string extension = fileName.substr(fileName.find_last_of('.') + 1);
        transform(extension.begin(), extension.end(), extension.begin(), ::tolower);
        // ignore temp files
        if (extension.compare("tmp") == 0) {
            continue;
        }

        if (strcmp(ptr->d_name, ".") == 0 || strcmp(ptr->d_name, "..") == 0){
            continue;
        }

        stat(path.c_str(), &s_buff);

        if (S_ISREG(s_buff.st_mode)) {
            filenames.push_back(fileName);
            remove(path.c_str());
        }

    }
    closedir(pDir);
}

void to_npy(std::string& outputFilename, float* array, unsigned int size, unsigned int channels) {
    const std::vector<unsigned long int> leshape11{3, size, size};
    std::vector<float> deit_vec(array, array + size * size * channels);
    npy::npy_data<float> data11;
    data11.data = deit_vec;
    data11.shape = leshape11;
    data11.fortran_order = false;
    std::string tmpFileName = outputFilename + ".tmp";
    write_npy(tmpFileName, data11);
    rename(tmpFileName.c_str(), outputFilename.c_str());
}
