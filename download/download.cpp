#include <iostream>
#include <cstring>
#include <curl/curl.h>
#include "Task.h"
#include "ThreadPool.h"
#include <random>
#include <algorithm>
#include "unistd.h"
#include "libcuda.h"

using namespace std;

static size_t write_data(void *ptr, size_t size, size_t nmemb, void *stream) {
    size_t written = fwrite(ptr, size, nmemb, (FILE*) stream);
    return written;
}

CURL* build_curl_handle() {
    CURL *curl_handle;

    curl_handle = curl_easy_init();

//    curl_easy_setopt(curl_handle, CURLOPT_VERBOSE, 1L);

    curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 1L);

    curl_easy_setopt(curl_handle, CURLOPT_FOLLOWLOCATION, 1L);

    curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYPEER, 0);

    curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYHOST, 0);

    curl_easy_setopt(curl_handle, CURLOPT_CAINFO, "ca-bundle.crt");

    curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);

    return curl_handle;
}

void downloadFile(const string& filename, const string& newFilename) {
    CURL *curl_handle;
    const char *pagefilename = newFilename.c_str();
    FILE *pagefile;
    char *p = (char *)filename.c_str();

//    curl_handle = (CURL*)safeMap.putIfabsent(build_curl_handle);
    curl_handle = build_curl_handle();

    curl_easy_setopt(curl_handle, CURLOPT_URL, p);

    pagefile = fopen(pagefilename, "wb");
    if (pagefile) {
        curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, pagefile);

        CURLcode code = curl_easy_perform(curl_handle);
        if (code != CURLE_OK)
            std::cout<<"tid: "<<std::this_thread::get_id()<<" failed to download image "<< filename <<std::endl;
        else
            std::cout<<"tid: "<<std::this_thread::get_id()<<" downloaded " << pagefilename <<std::endl;

        fclose(pagefile);
    } else {
        cout<<"failed to download open output file: <"<<pagefilename<<">"<<endl;
    }

    curl_easy_cleanup(curl_handle);

    return;
}


std::random_device dev;
std::mt19937 rng(dev());
std::uniform_int_distribution<std::mt19937::result_type> dist15(1,15); // distribution in range [1, 6]

std::string temp_dir;

void download_proc(std::string& image_suffix, SafeQueue* queue) {

    string output_file_path(temp_dir);
    output_file_path.append("/").append(image_suffix);

    string url("https://image");
    unsigned long r = dist15(rng);
    url.append(std::to_string(r));
    url.append(".coupangcdn.com/image/");

    std::string url_image_suffix = image_suffix;
    replace(url_image_suffix.begin(), url_image_suffix.end(), '_', '/');
    url.append(url_image_suffix);

    downloadFile(url, output_file_path);
    ImageBase* image_base = new ImageBase();
    image_base->image_path = output_file_path;
    image_base->image_suffix = image_suffix;

    int ret = load_image_file(image_base->image_path.c_str(), image_base);
    if (ret == -1) {
        onInvalidImage(&image_base->pixel_data_d, image_base->source_xsize, image_base->source_ysize, image_base->channels);
        std::cout<<"create empty image for " << image_base->image_path <<" as stub"<<std::endl;
    }

    queue->enqueue(image_base);
}


SafeQueue* result = nullptr;
ThreadPool* pool;

void init(const int thread_count, const string& _temp_dir) {
    curl_global_init(CURL_GLOBAL_ALL);
    pool = new ThreadPool(thread_count);
    pool->init();
    temp_dir = _temp_dir;
}

void submit_download_job(std::vector<std::string>& image_suffix_fixed_vec){
    if (result != nullptr) {
        delete result;
        result = nullptr;
    }

    result = new SafeQueue();

    for (vector<string>::iterator iter = image_suffix_fixed_vec.begin(); iter != image_suffix_fixed_vec.end(); iter++) {
        string image_suffix_fixed = *iter;
        replace(image_suffix_fixed.begin(), image_suffix_fixed.end(), '_', '/');
        Task t(*iter, download_proc, result);
        pool->submit(t);
    }
}

void wait_downloading_over(std::vector<std::string>& image_suffix_vec) {
    while (result->size() < image_suffix_vec.size()) {
        usleep(10 * 1000);
        std::cout<<"waiting"<<std::endl;
    }
}

void clean() {
    pool->shutdown();
    delete pool;
}

ImageBase* get_next_image_base() {
    if (result->size() == 0)
        return nullptr;

    ImageBase* image_base = result->dequeue();
    return image_base;
}