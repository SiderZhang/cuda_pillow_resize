//
// Created by siderzhangPC on 2024/4/21.
//

#ifndef LIBJPEGTEST_JPEG_H
#define LIBJPEGTEST_JPEG_H

int read_jpeg_file(const char *input_filename, unsigned char **output_buffer, unsigned int *width, unsigned int *height, unsigned int *channels);

#endif //LIBJPEGTEST_JPEG_H
