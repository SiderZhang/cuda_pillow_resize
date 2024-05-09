#include <iostream>
#include <cmath>
#include <stdio.h>
#include "jpeg.h"
#include <vector>
#include <string.h>
#include <dirent.h>
#include "npy.h"
#include "sys/types.h"
#include "sys/stat.h"
#include <algorithm>

#include <sys/timeb.h>
#include <unistd.h>

#define PRECISION_BITS (32 - 8 - 2)


typedef signed char         INT8, *PINT8;
typedef signed short        INT16, *PINT16;
typedef signed int          INT32, *PINT32;
typedef unsigned char       UINT8, *PUINT8;
typedef unsigned short      UINT16, *PUINT16;
typedef unsigned int        UINT32, *PUINT32;

__global__ void coeffs(int inSize, int in0, int in1,
                       int outSize, int *boundsp, double *kkp) {
    double scale, filterscale;
    double center, ww, ss;
    int x, ksize, xmin, xmax;
    int *bounds;
    double *kk, *k;

    /* prepare for horizontal stretch */
    filterscale = scale = (double)(in1 - in0) / outSize;
    if (filterscale < 1.0) {
        filterscale = 1.0;
    }

    /* determine support size (length of resampling filter) */
    double support = 1.0 * filterscale;

    /* maximum number of coeffs */
    ksize = (int)ceil(support) * 2 + 1;

    // check for overflow
    if (outSize > INT_MAX / (ksize * (int)sizeof(double))) {
        return;
    }

    /* coefficient buffer */
    /* malloc check ok, overflow checked above */
    kk = kkp;

    /* malloc check ok, ksize*sizeof(double) > 2*sizeof(int) */
    bounds = boundsp;
    int _xx = blockIdx.x * blockDim.x + threadIdx.x;
    if (_xx >= outSize) {
        return;
    }

    for (unsigned int xx = _xx;xx < outSize;xx += gridDim.x * blockDim.x) {
        center = in0 + (xx + 0.5) * scale;
        ww = 0.0;
        ss = 1.0 / filterscale;
        // Round the value
        xmin = (int) (center - support + 0.5);
        if (xmin < 0) {
            xmin = 0;
        }
        // Round the value
        xmax = (int) (center + support + 0.5);
        if (xmax > inSize) {
            xmax = inSize;
        }
        xmax -= xmin;
        k = &kk[xx * ksize];
        for (x = 0; x < xmax; x++) {
            double filteringX = (x + xmin - center + 0.5) * ss;

            double w = 0;//filterp->filter((x + xmin - center + 0.5) * ss);

            if (filteringX < 0.0) {
                filteringX = -filteringX;
            }
            if (filteringX < 1.0) {
                w = 1.0 - filteringX;
            } else {
                w = 0.0;
            }

            k[x] = w;
            ww += w;
        }
        for (x = 0; x < xmax; x++) {
            if (ww != 0.0) {
                k[x] /= ww;
            }
        }
        // Remaining values should stay empty if they are used despite of xmax.
        for (; x < ksize; x++) {
            k[x] = 0;
        }

        bounds[xx * 2 + 0] = xmin;
        bounds[xx * 2 + 1] = xmax;
    }
}


__global__ void normalize_coeffs(int outSize, int ksize, double *prekk) {
    INT32 *kk;

    // use the same buffer for normalized coefficients
    kk = (INT32 *)prekk;

    unsigned int _x = blockIdx.x * blockDim.x + threadIdx.x;
    if (_x >= outSize * ksize) {
        return;
    }


    for (unsigned int x = _x;x < outSize * ksize;x += gridDim.x * blockDim.x) {
        if (prekk[x] < 0) {
            kk[x] = (int) (-0.5 + prekk[x] * (1 << PRECISION_BITS));
        } else {
            kk[x] = (int) (0.5 + prekk[x] * (1 << PRECISION_BITS));
        }
    }
}

__global__ void shift_ysize(int *boundsp, int ysize) {
    unsigned int _i = blockIdx.x * blockDim.x + threadIdx.x;
    if (_i >= ysize) {
        return;
    }
//    printf("Bound %d %d\n", boundsp[i * 2], boundsp[0]);

    // Shift bounds for vertical pass
    for (unsigned int i = _i;i < ysize;i += gridDim.x * blockDim.x) {
        boundsp[i * 2] -= boundsp[0];
    }
}

__global__ void build_result_horiz(int srcXsize, int channelCount, int xsize, int ysize,
                                   unsigned char* input, unsigned char* output, int ksize, int *bounds, double *prekk, unsigned char* _lookups) {
    unsigned int _id = blockIdx.x * blockDim.x + threadIdx.x;
    if (_id >= xsize * ysize) {
        return;
    }

    for (unsigned int id = _id;id < xsize * ysize;id += gridDim.x * blockDim.x) {

        unsigned int xx = id % xsize;
        unsigned int yy = id / xsize;

        int ss0, ss1, ss2;
        INT32 *kk = (INT32 *) prekk;

        int xmin = bounds[xx * 2 + 0];
        int xmax = bounds[xx * 2 + 1];
        INT32 *k = &kk[xx * ksize];
        ss0 = ss1 = ss2 = 1 << (PRECISION_BITS - 1);

        for (int x = 0; x < xmax; x++) {

            ss0 += ((UINT8) input[yy * srcXsize * channelCount + (x + xmin) * channelCount + 0]) *
                   k[x];
            ss1 += ((UINT8) input[yy * srcXsize * channelCount + (x + xmin) * channelCount + 1]) *
                   k[x];
            ss2 += ((UINT8) input[yy * srcXsize * channelCount + (x + xmin) * channelCount + 2]) *
                   k[x];
        }

        UINT8 *lookups = &_lookups[640];
        UINT8 ss0_1 = lookups[ss0 >> PRECISION_BITS];
        UINT8 ss1_1 = lookups[ss1 >> PRECISION_BITS];
        UINT8 ss2_1 = lookups[ss2 >> PRECISION_BITS];

        output[yy * xsize * channelCount + xx * channelCount + 0] = ss0_1;
        output[yy * xsize * channelCount + xx * channelCount + 1] = ss1_1;
        output[yy * xsize * channelCount + xx * channelCount + 2] = ss2_1;
    }
}

__global__ void build_result_vert(int srcXsize, int channelCount, int xsize, int ysize,
                                  unsigned char* input, unsigned char* output, int ksize, int *bounds, double *prekk, unsigned char* _lookups) {
    unsigned int _id = blockIdx.x * blockDim.x + threadIdx.x;
    if (_id >= xsize * ysize) {
        return;
    }

    for (unsigned int id = _id;id < xsize * ysize;id += gridDim.x * blockDim.x) {
        unsigned int xx = id % xsize;
        unsigned int yy = id / xsize;

        int ss0, ss1, ss2;
        INT32 *kk = (INT32 *) prekk;

        INT32 *k = &kk[yy * ksize];
        int ymin = bounds[yy * 2 + 0];
        int ymax = bounds[yy * 2 + 1];

        ss0 = ss1 = ss2 = 1 << (PRECISION_BITS - 1);
        for (int y = 0; y < ymax; y++) {
            ss0 += ((UINT8) input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 0]) * k[y];
            ss1 += ((UINT8) input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 1]) * k[y];
            ss2 += ((UINT8) input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 2]) * k[y];
        }

        UINT8 *lookups = &_lookups[640];
        UINT8 ss0_1 = lookups[ss0 >> PRECISION_BITS];
        UINT8 ss1_1 = lookups[ss1 >> PRECISION_BITS];
        UINT8 ss2_1 = lookups[ss2 >> PRECISION_BITS];
        output[yy * xsize * channelCount + xx * channelCount + 0] = ss0_1;
        output[yy * xsize * channelCount + xx * channelCount + 1] = ss1_1;
        output[yy * xsize * channelCount + xx * channelCount + 2] = ss2_1;
    }
}

int *bounds_horiz_d;
double *kk_horiz_d;
int ksize_horiz;

int  *bounds_vert_d;
double *kk_vert_d;
int ksize_vert;

void coffes(unsigned int im_xsize, unsigned int im_ysize, int xsize, int ysize) {

    float box[4] = {0, 0, 1.0f * im_xsize, 1.0f * im_ysize};

    double filterscale_horiz = (double)(box[2] - box[0]) / xsize;
    if (filterscale_horiz < 1.0) {
        filterscale_horiz = 1.0;
    }
    ksize_horiz = (int)ceil(filterscale_horiz) * 2 + 1;

    cudaMalloc(&kk_horiz_d, xsize * ksize_horiz * sizeof(double));
    cudaMalloc(&bounds_horiz_d, xsize * 2 * sizeof(int));

    coeffs<<<256, 256>>>(im_xsize,
                        box[0],
                        box[2],
                        xsize,
                        bounds_horiz_d,
                        kk_horiz_d);

    normalize_coeffs<<<xsize, ksize_horiz>>>(xsize, ksize_horiz, kk_horiz_d);

    double filterscale_vert = (double)(box[3] - box[1]) / ysize;
    if (filterscale_vert < 1.0) {
        filterscale_vert = 1.0;
    }
    ksize_vert = (int)ceil(filterscale_vert) * 2 + 1;
    cudaMalloc(&kk_vert_d, ysize * ksize_vert * sizeof(double));
    cudaMalloc(&bounds_vert_d, ysize * 2 * sizeof(int));

    coeffs<<<256, 256>>>(im_ysize,
                        box[1],
                        box[3],
                        ysize,
                        bounds_vert_d,
                        kk_vert_d);

    normalize_coeffs<<<ysize, ksize_vert>>>(ysize, ksize_vert, kk_vert_d);

    shift_ysize<<<256, 256>>>(bounds_vert_d, ysize);
}

__global__ void rescale_normalize_d(unsigned char *data, float* result, int xsize, int ysize) {
    unsigned int _id = blockIdx.x * blockDim.x + threadIdx.x;

    if (_id >= xsize * ysize) {
        return;
    }

    for (unsigned int id = _id;id < xsize * ysize;id += gridDim.x * blockDim.x) {
        // rescale
        double pixel0 = ((double)data[id * 3 + 0]) * 1.0 / 255;
        double pixel1 = ((double)data[id * 3 + 1]) * 1.0 / 255;
        double pixel2 = ((double)data[id * 3 + 2]) * 1.0 / 255;

        // normalize
        pixel0 = (pixel0 - 0.5) / 0.5;
        pixel1 = (pixel1 - 0.5) / 0.5;
        pixel2 = (pixel2 - 0.5) / 0.5;

        // transpose
        result[id + 0 * xsize * ysize] = (float)pixel0;
        result[id + 1 * xsize * ysize] = (float)pixel1;
        result[id + 2 * xsize * ysize] = (float)pixel2;
    }
}

UINT8 lookups_h[1280] = {
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
        0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   1,   2,   3,   4,   5,
        6,   7,   8,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,
        23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,
        40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,
        57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,
        74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
        91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107,
        108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124,
        125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141,
        142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158,
        159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
        176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192,
        193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
        210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226,
        227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243,
        244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255,
};

extern "C" int resize(const char* filename, std::string& outputFilename, unsigned int target_xsize, unsigned int target_ysize) {
    unsigned char *image_pixel_data = NULL;
    unsigned int source_xsize;
    unsigned int source_ysize;
    unsigned int channels;

    int ret = read_jpeg_file(filename, &image_pixel_data, &source_xsize, &source_ysize, &channels);
    if (ret != 0) {
        std::cerr<<"failed to load file" << filename <<std::endl;
        return -1;
    }

    coffes(source_xsize, source_ysize, target_xsize, target_ysize);

    unsigned char* input;
    cudaMalloc(&input, source_xsize * source_ysize * channels * sizeof(unsigned char));

    cudaMemcpy(input, image_pixel_data, source_xsize * source_ysize * channels * sizeof(unsigned char), cudaMemcpyHostToDevice);

    unsigned char* temp;
    cudaMalloc(&temp, target_xsize * source_ysize * channels * sizeof(unsigned char));

    unsigned char* result_pixel_data_d;
    cudaMalloc(&result_pixel_data_d, target_xsize * target_ysize * channels * sizeof(unsigned char));

    unsigned char* looksup_d;
    cudaMalloc(&looksup_d, 1280 * sizeof(unsigned char));
    cudaMemcpy(looksup_d, lookups_h, 1280 * sizeof(unsigned char), cudaMemcpyHostToDevice);

    build_result_horiz<<<256, 1024>>>(source_xsize, channels, target_xsize, source_ysize, input, temp, ksize_horiz, bounds_horiz_d, kk_horiz_d, looksup_d);
    build_result_vert<<<256, 1024>>>(target_xsize, channels, target_xsize, target_ysize, temp, result_pixel_data_d, ksize_vert, bounds_vert_d, kk_vert_d, looksup_d);

    float* normalized_pixel_data_d;
    cudaMalloc(&normalized_pixel_data_d, target_xsize * target_ysize * channels * sizeof(float));

    rescale_normalize_d<<<256, 1024>>>(result_pixel_data_d, normalized_pixel_data_d, target_xsize, target_ysize);

    std::vector<float> v;
    v.resize(target_xsize * target_ysize * 3);
    float* normalized_pixel_data = (float*) calloc(target_xsize * target_ysize, channels * sizeof(float));
    cudaMemcpy(&v[0], normalized_pixel_data_d, target_xsize * target_ysize * channels * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(bounds_horiz_d);
    cudaFree(kk_horiz_d);
    cudaFree(bounds_vert_d);
    cudaFree(kk_vert_d);

    cudaFree(input);
    cudaFree(temp);
    cudaFree(looksup_d);
    cudaFree(result_pixel_data_d);
    cudaFree(normalized_pixel_data_d);

    cudaDeviceSynchronize();


    const std::vector<unsigned long> leshape11{3, target_xsize, target_ysize};

    const npy::npy_data<float> data11{v, leshape11, false};
    std::string tmpFileName = outputFilename + ".tmp";
    write_npy(tmpFileName, data11);
    rename(tmpFileName.c_str(), outputFilename.c_str());

    return 0;
}

void readDir(const char* dirPath, std::vector<std::string>& filenames) {
    DIR *pDir;
    struct dirent* ptr;
    if(!(pDir = opendir(dirPath))){
        std::cout<<"Folder doesn't Exist!"<<std::endl;
        return;
    }
    struct stat s_buff;

    while((ptr = readdir(pDir)) != nullptr) {
        std::string fileName(ptr->d_name);
        std::string path = std::string(dirPath) + "/" + fileName;

        std::string extension = fileName.substr(fileName.find_last_of(".") + 1);
        transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

        if (!extension.compare("jpg") == 0 && !extension.compare("jpeg") == 0) {
            continue;
        }

        if (strcmp(ptr->d_name, ".") == 0 || strcmp(ptr->d_name, "..") == 0){
            continue;
        }

        stat(path.c_str(), &s_buff);
        if (!S_ISREG(s_buff.st_mode)) {
            continue;
        }

        filenames.push_back(fileName);
    }
    closedir(pDir);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        std::cerr<<"arguments: dirPath width height outputfilename"<<std::endl;
        return -1;
    }
    char* dirPath = argv[1];
    char* outputFilename = argv[2];
    unsigned int xsize = std::atoi(argv[3]);
    unsigned int ysize = std::atoi(argv[4]);

    while (true) {
        std::vector<std::string> fileNames;
        readDir(dirPath, fileNames);

        timeb t;
        ftime(&t);
        long t1 = t.time * 1000 + t.millitm;
        int count = 0;
        for (auto iter = fileNames.begin();iter != fileNames.end();iter++) {
            std::string inputFilename22 = std::string(dirPath) + "/" + *iter;
            std::string outputFilename22 = std::string(outputFilename) + "/" + *iter + ".npy";

            FILE *output_file_test;
            if ((output_file_test = fopen(outputFilename22.c_str(), "rb")) != NULL) {
                std::cerr<<"npy file exists <" << outputFilename22 << ">" << " skip it" <<std::endl;
                fclose(output_file_test);
                continue;
            }

            int ret = resize(inputFilename22.c_str(), outputFilename22, xsize, ysize);
            if (ret != 0)
                continue;

            count++;
        }
        ftime(&t);
        long t2 = t.time * 1000 + t.millitm;

        if (count != 0) {
            std::cout << "process iamges " << count << " for time " << t2 - t1 << " millis at "<< t2 << std::endl;
        }
        usleep(5000);
    }
    return 0;
}
