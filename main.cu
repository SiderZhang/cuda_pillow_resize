#include <iostream>
#include <cmath>
#include <stdio.h>
#include <vector>
#include <string.h>
#include "npy.h"

#include "libcuda.h"

#include <sys/timeb.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/stat.h>

#define PRECISION_BITS (32 - 8 - 2)


typedef signed char         INT8, *PINT8;
typedef signed short        INT16, *PINT16;
typedef signed int          INT32, *PINT32;
typedef unsigned char       UINT8, *PUINT8;
typedef unsigned short      UINT16, *PUINT16;
typedef unsigned int        UINT32, *PUINT32;

#define LANZOS_INTERPOLATION 1
#define BILINEAR_INTERPOLATION 2
#define BICUBIC_INTERPOLATION 3
#define BOX_INTERPOLATION 4
#define HAMMING_INTERPOLATION 5

struct filter {
    double (*filter)(double x);
    double support;
};

__device__ __forceinline__ double
box_filter(double x) {
    if (x > -0.5 && x <= 0.5) {
        return 1.0;
    }
    return 0.0;
}

__device__ double
bilinear_filter(double x) {
    if (x < 0.0) {
        x = -x;
    }
    if (x < 1.0) {
        return 1.0 - x;
    }
    return 0.0;
}

__device__ double
hamming_filter(double x) {
    if (x < 0.0) {
        x = -x;
    }
    if (x == 0.0) {
        return 1.0;
    }
    if (x >= 1.0) {
        return 0.0;
    }
    x = x * M_PI;
    return sin(x) / x * (0.54f + 0.46f * cos(x));
}

__device__ double
bicubic_filter(double x) {
    /* https://en.wikipedia.org/wiki/Bicubic_interpolation#Bicubic_convolution_algorithm
     */
#define a -0.5
    if (x < 0.0) {
        x = -x;
    }
    if (x < 1.0) {
        return ((a + 2.0) * x - (a + 3.0)) * x * x + 1;
    }
    if (x < 2.0) {
        return (((x - 5) * x + 8) * x - 4) * a;
    }
    return 0.0;
#undef a
}

__device__ __forceinline__ double
sinc_filter(double x) {
    if (x == 0.0) {
        return 1.0;
    }
    x = x * M_PI;
    return sin(x) / x;
}

__device__ __forceinline__ double
lanczos_filter(double x) {
    /* truncated sinc */
    if (-3.0 <= x && x < 3.0) {
        return sinc_filter(x) * sinc_filter(x / 3);
    }
    return 0.0;
}


__device__ static struct filter BOX_D = {box_filter, 0.5};
__device__ static struct filter BILINEAR_D = {bilinear_filter, 1.0};
__device__ static struct filter HAMMING_D = {hamming_filter, 1.0};
__device__ static struct filter BICUBIC_D = {bicubic_filter, 2.0};
__device__ static struct filter LANCZOS_D = {lanczos_filter, 3.0};

static struct filter BOX_H = {box_filter, 0.5};
static struct filter BILINEAR_H = {bilinear_filter, 1.0};
static struct filter HAMMING_H = {hamming_filter, 1.0};
static struct filter BICUBIC_H = {bicubic_filter, 2.0};
static struct filter LANCZOS_H = {lanczos_filter, 3.0};

__global__ void test(int step, int* boundsp, int size) {
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= size) {
        return;
    }

    printf("step: %d, id: %d, bound: %d\n", step, id, boundsp[id]);
}


__global__ void coeffs(int inSize, int in0, int in1,
                       int outSize, int *boundsp, double *kkp, int interpolation_mode) {
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
    struct filter *filterp = NULL;

    switch (interpolation_mode){
        case 1:
            filterp = &LANCZOS_D;
            break;
        case 2:
            filterp = &BILINEAR_D;
            break;
        case 3:
            filterp = &BICUBIC_D;
            break;
        case 4:
            filterp = &BOX_D;
            break;
        case 5:
            filterp = &HAMMING_D;
            break;
        default:
            printf("unsupported interpolation mode %d", interpolation_mode);
            return;
    }

    double support = filterp->support * filterscale;
//    double support = supportFacter * filterscale;

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

//            double filteringX = (x + xmin - center + 0.5) * ss;
//
//            double w = 0;//filterp->filter((x + xmin - center + 0.5) * ss);
//
//            if (filteringX < 0.0) {
//                filteringX = -filteringX;
//            }
//            if (filteringX < 1.0) {
//                w = 1.0 - filteringX;
//            } else {
//                w = 0.0;
//            }
            double w = filterp->filter((x + xmin - center + 0.5) * ss);
//            double w = pFilterFunc((x + xmin - center + 0.5) * ss);

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


__global__ void normalize_coeffs(unsigned int outSize, int ksize, double *prekk) {
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

__global__ void build_result_horiz(unsigned int srcXsize, unsigned int channelCount, unsigned int xsize, unsigned int ysize,
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
//            int s = yy * srcXsize * channelCount + (x + xmin) * channelCount;
//            if (s < 0 || s >= srcXsize * ysize * channelCount) {
//                printf("invalid index: %d, %d, %d, %d, %d, %d\n", yy, channelCount, x, xmax, xmin, s);
//            }

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

__global__ void build_result_vert(unsigned int srcXsize, unsigned int channelCount, unsigned int xsize, unsigned int ysize,
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

float* corp(float* src, unsigned int input_size, unsigned int corp_size, unsigned int channels){
    float* target;
    cudaMalloc(&target, corp_size * corp_size * channels * sizeof(float));
    unsigned int left = (input_size - corp_size) / 2;
    unsigned int top = (input_size - corp_size) / 2;
    for (int i = 0;i < corp_size;i++) {
        unsigned int src_start = left + input_size * (i + top);
        unsigned int dest_start = corp_size * i;
        cudaMemcpy(&target[dest_start + corp_size * corp_size * 0], &src[src_start + input_size * input_size * 0], corp_size * sizeof(float), cudaMemcpyDeviceToDevice);
        cudaMemcpy(&target[dest_start + corp_size * corp_size * 1], &src[src_start + input_size * input_size * 1], corp_size * sizeof(float), cudaMemcpyDeviceToDevice);
        cudaMemcpy(&target[dest_start + corp_size * corp_size * 2], &src[src_start + input_size * input_size * 2], corp_size * sizeof(float), cudaMemcpyDeviceToDevice);
    }
    cudaFree(src);
    return target;
}

__global__ void rescale_normalize_d(unsigned char *data, float* result, unsigned int xsize, unsigned int ysize,
                                    float mean0, float mean1, float mean2, float std0, float std1, float std2) {
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
        pixel0 = (pixel0 - mean0) / std0;
        pixel1 = (pixel1 - mean1) / std1;
        pixel2 = (pixel2 - mean2) / std2;

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
float* image_preprocess(unsigned char* input, unsigned int source_xsize, unsigned int source_ysize, unsigned int target_xsize, unsigned int target_ysize, unsigned int channels, int interpolation_mode, float* mean, float* std) {
    int *bounds_horiz_d;
    double *kk_horiz_d;
    int ksize_horiz;

    int  *bounds_vert_d;
    double *kk_vert_d;
    int ksize_vert;

    float box[4] = {0, 0, 1.0f * source_xsize, 1.0f * source_ysize};

    double support;
    struct filter *filterp;

    switch (interpolation_mode){
        case 1:
            filterp = &LANCZOS_H;
            break;
        case 2:
            filterp = &BILINEAR_H;
            break;
        case 3:
            filterp = &BICUBIC_H;
            break;
        case 4:
            filterp = &BOX_H;
            break;
        case 5:
            filterp = &HAMMING_H;
            break;
        default:
            printf("unsupported interpolation mode %d", interpolation_mode);
            return NULL;
    }

    double filterscale_horiz = (double)(box[2] - box[0]) / target_xsize;
    if (filterscale_horiz < 1.0) {
        filterscale_horiz = 1.0;
    }
    ksize_horiz = (int)ceil(filterscale_horiz * filterp->support) * 2 + 1;

    cudaMalloc(&kk_horiz_d, target_xsize * ksize_horiz * sizeof(double));
    cudaMalloc(&bounds_horiz_d, target_xsize * 2 * sizeof(int));

    coeffs<<<256, 256>>>(source_xsize,
                         box[0],
                         box[2],
                         target_xsize,
                         bounds_horiz_d,
                         kk_horiz_d,
                         interpolation_mode);

    cudaError_t error =  cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    normalize_coeffs<<<target_xsize, ksize_horiz>>>(target_xsize, ksize_horiz, kk_horiz_d);

    error =  cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    double filterscale_vert = (double)(box[3] - box[1]) / target_ysize;
    if (filterscale_vert < 1.0) {
        filterscale_vert = 1.0;
    }
    ksize_vert = (int)ceil(filterscale_vert * filterp->support) * 2 + 1;
    cudaMalloc(&kk_vert_d, target_ysize * ksize_vert * sizeof(double));
    cudaMalloc(&bounds_vert_d, target_ysize * 2 * sizeof(int));

    error = cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    coeffs<<<256, 256>>>(source_ysize,
                         box[1],
                         box[3],
                         target_ysize,
                         bounds_vert_d,
                         kk_vert_d,
                         interpolation_mode);

    error = cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    normalize_coeffs<<<target_ysize, ksize_vert>>>(target_ysize, ksize_vert, kk_vert_d);

    error = cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    shift_ysize<<<256, 256>>>(bounds_vert_d, target_ysize);

    error = cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    unsigned char* temp;
    cudaMalloc(&temp, target_xsize * source_ysize * channels * sizeof(unsigned char));

    unsigned char* result_pixel_data_d;
    cudaMalloc(&result_pixel_data_d, target_xsize * target_ysize * channels * sizeof(unsigned char));

    unsigned char* looksup_d;
    cudaMalloc(&looksup_d, 1280 * sizeof(unsigned char));
    cudaMemcpy(looksup_d, lookups_h, 1280 * sizeof(unsigned char), cudaMemcpyHostToDevice);

    build_result_horiz<<<256, 1024>>>(source_xsize, channels, target_xsize, source_ysize, input, temp, ksize_horiz, bounds_horiz_d, kk_horiz_d, looksup_d);

    error = cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    build_result_vert<<<256, 1024>>>(target_xsize, channels, target_xsize, target_ysize, temp, result_pixel_data_d, ksize_vert, bounds_vert_d, kk_vert_d, looksup_d);

    error = cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    float* normalized_pixel_data_d;
    cudaMalloc(&normalized_pixel_data_d, target_xsize * target_ysize * channels * sizeof(float));

    rescale_normalize_d<<<256, 1024>>>(result_pixel_data_d, normalized_pixel_data_d, target_xsize, target_ysize,
                                       mean[0], mean[1], mean[2], std[0], std[1], std[2]);

    error = cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    cudaFree(temp);
    cudaFree(looksup_d);
    cudaFree(result_pixel_data_d);

    cudaFree(bounds_horiz_d);
    cudaFree(kk_horiz_d);
    cudaFree(bounds_vert_d);
    cudaFree(kk_vert_d);

    error = cudaGetLastError();
    if (error != cudaError_t::cudaSuccess) {
        std::cerr<<cudaGetErrorString(error)<<std::endl;
        return NULL;
    }

    return normalized_pixel_data_d;
}

float* vit_preprocess(unsigned char* input, unsigned int source_xsize, unsigned int source_ysize, unsigned int vit_size, unsigned int channels) {
    float mean[3] = {0.5, 0.5, 0.5};
    float std[3] = {0.5, 0.5, 0.5};
    int interpolation_mode = BILINEAR_INTERPOLATION;
    return image_preprocess(input, source_xsize, source_ysize, vit_size, vit_size, channels, interpolation_mode, mean, std);
}

float* deit_preprocess(unsigned char* input, unsigned int source_xsize, unsigned int source_ysize, unsigned int vit_size, unsigned int corp_size, unsigned int channels) {
    float mean[3] = {0.485, 0.456, 0.406};
    float std[3] = {0.229, 0.224, 0.225};
    int interpolation_mode = BICUBIC_INTERPOLATION;
    float* resizedImage = image_preprocess(input, source_xsize, source_ysize, vit_size, vit_size, channels, interpolation_mode, mean, std);

    return resizedImage;
//    return corp(resizedImage, vit_size, corp_size, channels);
}

float* vis_preprocess(unsigned char* input, unsigned int source_xsize, unsigned int source_ysize, unsigned int vit_size, unsigned int channels) {
    float mean[3] = {0.48145466, 0.4578275, 0.40821073};
    float std[3] = {0.26862954, 0.26130258, 0.27577711};
    int interpolation_mode = BICUBIC_INTERPOLATION;
    float* resizedImage = image_preprocess(input, source_xsize, source_ysize, vit_size, vit_size, channels, interpolation_mode, mean, std);
    return resizedImage;
}

void onImageRead(unsigned char* input, unsigned char **output_buffer, unsigned int width, unsigned int height, unsigned int channels){
    cudaMalloc(output_buffer, width * height * channels * sizeof(unsigned char));
    cudaMemcpy(*output_buffer, input, width * height * channels * sizeof(unsigned char), cudaMemcpyHostToDevice);
}

void to_npy(const char* filename, std::string& outputFilename, float* array, unsigned int size, unsigned int channels) {
    const std::vector<unsigned long int> leshape11{3, size, size};
    std::vector<float> deit_vec(array, array + size * size * channels);
//    std::vector<float> deit_vec;
    npy::npy_data<float> data11;
    data11.data = deit_vec;
    data11.shape = leshape11;
    data11.fortran_order = false;
    std::string tmpFileName = outputFilename + ".tmp";
    write_npy(tmpFileName, data11);
    rename(tmpFileName.c_str(), outputFilename.c_str());
}


int image_process(const char* filename, std::string& output_filename_prefix, unsigned int vit_size, unsigned int vis_size, unsigned int deit_size, unsigned int corp_size) {
    unsigned int source_xsize;
    unsigned int source_ysize;
    unsigned int channels;

    unsigned char* input;
    int ret = load_image_file(filename, &input, source_xsize, source_ysize, channels);
    if (ret == -1) {
        return -1;
    }

    float* deit_pixel_data_d = deit_preprocess(input, source_xsize, source_ysize, deit_size, corp_size, channels);
    float* deit_result_h;
    deit_result_h = (float*)malloc(corp_size * corp_size * channels * sizeof(float));
    cudaMemcpy(deit_result_h, deit_pixel_data_d, corp_size * corp_size * channels * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(deit_pixel_data_d);
    std::string deit_output_filename = output_filename_prefix + "_deit.npy";
    to_npy(filename, deit_output_filename, deit_result_h, corp_size, channels);
    free(deit_result_h);

    float* vis_pixel_data_d = vis_preprocess(input, source_xsize, source_ysize, vis_size, channels);
    float* vis_reuslt_h;
    vis_reuslt_h = (float*)malloc(vis_size * vis_size * channels * sizeof(float));
    cudaMemcpy(vis_reuslt_h, vis_pixel_data_d, vis_size * vis_size * channels * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(vis_pixel_data_d);
    std::string vis_output_filename = output_filename_prefix + "_vis.npy";
    to_npy(filename, vis_output_filename, vis_reuslt_h, vis_size, channels);
    free(vis_reuslt_h);

    float* vit_pixel_data_d = vit_preprocess(input, source_xsize, source_ysize, vit_size, channels);
    float* vit_reuslt_h;
    vit_reuslt_h = (float*)malloc(vit_size * vit_size * channels * sizeof(float));
    cudaMemcpy(vit_reuslt_h, vit_pixel_data_d, vit_size * vit_size * channels * sizeof(float), cudaMemcpyDeviceToHost);
    cudaFree(vit_pixel_data_d);
    std::string vit_output_filename = output_filename_prefix + "_vit.npy";
    to_npy(filename, vit_output_filename, vit_reuslt_h, vit_size, channels);
    free(vit_reuslt_h);

    cudaFree(input);
    cudaDeviceSynchronize();

    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr<<"arguments: input_dir output_dir"<<std::endl;
        return -1;
    }
    char* input_dir = argv[1];
    char* output_dir = argv[2];

    while (true) {
        std::vector<std::string> fileNames;
        readDir(input_dir, fileNames, std::string(output_dir));

        timeb t;
        ftime(&t);
        long t1 = t.time * 1000 + t.millitm;
        int count = 0;
        for (auto iter = fileNames.begin();iter != fileNames.end();iter++) {
            std::string input_file_abs_path = std::string(input_dir) + "/" + *iter;
            std::string output_file_abs_path = std::string(output_dir) + "/" + *iter;

            int ret = image_process(input_file_abs_path.c_str(), output_file_abs_path, 224, 224, 256, 256);
            if (ret != 0)
                continue;

            ret = remove(input_file_abs_path.c_str());
            if (ret != 0) {
                std::cerr<<"failed to delete processed image file "<<input_file_abs_path.c_str()<<std::endl;
                continue;
            }

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
