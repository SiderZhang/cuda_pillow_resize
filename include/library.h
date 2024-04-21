#ifndef LIBJPEGOP_LIBRARY_H
#define LIBJPEGOP_LIBRARY_H
int read_jpeg_file(const char *input_filename, unsigned char **output_buffer, unsigned int *width, unsigned int *height, unsigned int *channels);
#endif //LIBJPEGOP_LIBRARY_H
