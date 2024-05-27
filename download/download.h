//
// Created by siderzhangPC on 2024/5/23.
//

#ifndef CURLDOWNLOAD_DOWNLOAD_H
#define CURLDOWNLOAD_DOWNLOAD_H

#include <string>
#include <vector>

void init(const int thread_count, const std::string& _temp_dir);
void submit_download_job(std::vector<std::string>& image_suffix_fixed_vec);
void wait_downloading_over(std::vector<std::string>& image_suffix_vec);
void clean();
ImageBase* get_next_image_base();

#endif //CURLDOWNLOAD_DOWNLOAD_H