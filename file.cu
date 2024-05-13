#include <string>
#include <vector>
#include <iostream>

#include <algorithm>
#include <dirent.h>
#include "sys/stat.h"

using namespace std;

void scanDir(const string& base_path, const string& dir_path, std::vector<std::string>& filenames) {
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(std::string(base_path).append(dir_path).c_str()))){
        std::cout<<"Folder doesn't Exist!"<<std::endl;
        return;
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
            scanDir(base_path, relative_path, filenames);
            continue;
        }
    }
    closedir(pDir);
}

void readDir(const char* dirPath, std::vector<std::string>& filenames) {

    scanDir(std::string(dirPath), std::string(""), filenames);

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