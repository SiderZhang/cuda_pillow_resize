//
// Created by siderzhangPC on 2024/4/21.
//
#include <stdio.h>
#include <stdlib.h>
#include <jpeglib.h>
#include <cerrno>
#include <cstring>
#include <exception>
#include <iostream>
#include <stdexcept>

void jpeg_error(j_common_ptr cinfo) {
    char err_msg[JMSG_LENGTH_MAX];
    (*cinfo->err->format_message)(cinfo, err_msg);
    throw std::invalid_argument(err_msg);
}

int read_jpeg_file(const char *input_filename, unsigned char **output_buffer, unsigned int *width, unsigned int *height, unsigned int *channels) {
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;
    FILE *input_file;
    FILE *output_file;
    JSAMPARRAY buffer;
    int row_width;

    unsigned char *rowdata = NULL;

    if ((input_file = fopen(input_filename, "rb")) == NULL) {
        int err = errno;
        fprintf(stderr, "can't open %s, error: %s\n", input_filename, strerror(err));
        return -1;
    }
    cinfo.err = jpeg_std_error(&jerr);
    jerr.error_exit = jpeg_error;

    try {
        jpeg_create_decompress(&cinfo);
        jpeg_stdio_src(&cinfo, input_file);

        jpeg_read_header(&cinfo, TRUE);
        jpeg_start_decompress(&cinfo);
        row_width = cinfo.output_width * cinfo.output_components;
        *width = cinfo.output_width;
        *height = cinfo.output_height;
        *channels = cinfo.output_components;

        buffer = (*cinfo.mem->alloc_sarray)((j_common_ptr)&cinfo, JPOOL_IMAGE, row_width, 1);

        *output_buffer = (unsigned char*)malloc(row_width * cinfo.output_height);
        memset(*output_buffer, 0, row_width * cinfo.output_height);
        rowdata = *output_buffer;

        while (cinfo.output_scanline < cinfo.output_height) {
            jpeg_read_scanlines(&cinfo, buffer, 1);

            memcpy(rowdata, *buffer, row_width);
            rowdata += row_width;
        }

        jpeg_finish_decompress(&cinfo);
        jpeg_destroy_decompress(&cinfo);
    } catch(std::exception& e) {
        std::cerr<<e.what()<<std::endl;
        fclose(input_file);
        return -1;
    }
    fclose(input_file);

    return 0;
}