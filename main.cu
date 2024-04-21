#include <iostream>
#include <cmath>
#include "ipl.h"
#include <stdio.h>
#include "Imaging.h"
#include "include/library.h"

__global__ void hello() {
    printf("hello gpu\n");
}

#define MAKE_UINT32(u0, u1, u2, u3) \
    ((UINT32)(u0) | ((UINT32)(u1) << 8) | ((UINT32)(u2) << 16) | ((UINT32)(u3) << 24))

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
    int xx = blockIdx.x * blockDim.x + threadIdx.x;
    if (xx >= outSize) {
        return;
    }

        center = in0 + (xx + 0.5) * scale;
        ww = 0.0;
        ss = 1.0 / filterscale;
        // Round the value
        xmin = (int)(center - support + 0.5);
        if (xmin < 0) {
            xmin = 0;
        }
        // Round the value
        xmax = (int)(center + support + 0.5);
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


__global__ void normalize_coeffs(int outSize, int ksize, double *prekk) {
    INT32 *kk;

    // use the same buffer for normalized coefficients
    kk = (INT32 *)prekk;

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    if (x >= outSize * ksize) {
        return;
    }

    if (prekk[x] < 0) {
        kk[x] = (int)(-0.5 + prekk[x] * (1 << PRECISION_BITS));
    } else {
        kk[x] = (int)(0.5 + prekk[x] * (1 << PRECISION_BITS));
    }
}

__global__ void shift_ysize(int *boundsp, int ysize) {
    unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= ysize) {
        return;
    }
//    printf("Bound %d %d\n", boundsp[i * 2], boundsp[0]);

    // Shift bounds for vertical pass
    boundsp[i * 2] -= boundsp[0];
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

__global__ void build_result_horiz(int srcXsize, int channelCount, int xsize, int ysize,
                                   unsigned char* input, unsigned char* output, int ksize, int *bounds, double *prekk, unsigned char* _lookups) {
    unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= xsize * ysize) {
        return;
    }

    unsigned int xx = id % xsize;
    unsigned int yy = id / xsize;

    int ss0, ss1, ss2;
    INT32 *kk = (INT32 *)prekk;

    int xmin = bounds[xx * 2 + 0];
    int xmax = bounds[xx * 2 + 1];
    INT32 *k = &kk[xx * ksize];
    ss0 = ss1 = ss2 = 1 << (PRECISION_BITS - 1);

    for (int x = 0; x < xmax; x++) {

        ss0 += ((UINT8)input[yy * srcXsize * channelCount +  (x + xmin) * channelCount + 0]) *
               k[x];
        ss1 += ((UINT8)input[yy * srcXsize * channelCount +  (x + xmin) * channelCount + 1]) *
               k[x];
        ss2 += ((UINT8)input[yy * srcXsize * channelCount +  (x + xmin) * channelCount + 2]) *
               k[x];
//        if (yy == 0 && xx == 20) {
//            printf("%d %d\n index: %d %d %d\nsource: %d %d %d\n%d\n",
//                   xx, x,
//                   yy * srcXsize * channelCount +  (x + xmin) * channelCount + 0,
//                   yy * srcXsize * channelCount +  (x + xmin) * channelCount + 1,
//                   yy * srcXsize * channelCount +  (x + xmin) * channelCount + 2,
//                   input[yy * srcXsize * channelCount +  (x + xmin) * channelCount + 0],
//                   input[yy * srcXsize * channelCount +  (x + xmin) * channelCount + 1],
//                   input[yy * srcXsize * channelCount +  (x + xmin) * channelCount + 2],
//                   k[x]
//            );
//        }
    }

    UINT8 *lookups = &_lookups[640];
    UINT8 ss0_1 = lookups[ss0 >> PRECISION_BITS];
    UINT8 ss1_1 = lookups[ss1 >> PRECISION_BITS];
    UINT8 ss2_1 = lookups[ss2 >> PRECISION_BITS];
//    if (yy == 0 && xx == 20) {
//        printf("source ss: %d %d %d\nsource dd: %d %d %d\n",
//                ss0 >> PRECISION_BITS,
//                ss1 >> PRECISION_BITS,
//                ss2 >> PRECISION_BITS,
//                ss0_1, ss1_1, ss2_1);
//    }

    output[yy * xsize * channelCount + xx * channelCount + 0] = ss0_1;
    output[yy * xsize * channelCount + xx * channelCount + 1] = ss1_1;
    output[yy * xsize * channelCount + xx * channelCount + 2] = ss2_1;
}

__global__ void build_result_vert(int srcXsize, int channelCount, int xsize, int ysize,
                                  unsigned char* input, unsigned char* output, int ksize, int *bounds, double *prekk, unsigned char* _lookups) {
    unsigned int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id >= xsize * ysize) {
        return;
    }

    unsigned int xx = id % xsize;
    unsigned int yy = id / xsize;

    int ss0, ss1, ss2;
    INT32 *kk = (INT32 *)prekk;

    INT32 *k = &kk[yy * ksize];
    int ymin = bounds[yy * 2 + 0];
    int ymax = bounds[yy * 2 + 1];

    ss0 = ss1 = ss2 = 1 << (PRECISION_BITS - 1);
    for (int y = 0; y < ymax; y++) {
        ss0 += ((UINT8)input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 0]) * k[y];
        ss1 += ((UINT8)input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 1]) * k[y];
        ss2 += ((UINT8)input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 2]) * k[y];
        if (yy == 0 && xx == 0) {
            printf("%d %d\n index: %d %d %d\nsource: %d %d %d\n%d\n",
                   xx, y,
                   (y + ymin) * srcXsize * channelCount + xx * channelCount + 0,
                   (y + ymin) * srcXsize * channelCount + xx * channelCount + 1,
                   (y + ymin) * srcXsize * channelCount + xx * channelCount + 2,
                   input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 0],
                   input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 1],
                   input[(y + ymin) * srcXsize * channelCount + xx * channelCount + 2],
                   k[y]
            );
        }
    }

    UINT8 *lookups = &_lookups[640];
    UINT8 ss0_1 = lookups[ss0 >> PRECISION_BITS];
    UINT8 ss1_1 = lookups[ss1 >> PRECISION_BITS];
    UINT8 ss2_1 = lookups[ss2 >> PRECISION_BITS];
    output[yy * xsize * channelCount + xx * channelCount + 0] = ss0_1;
    output[yy * xsize * channelCount + xx * channelCount + 1] = ss1_1;
    output[yy * xsize * channelCount + xx * channelCount + 2] = ss2_1;
}

int *bounds_horiz_h, *bounds_horiz_d;
double *kk_horiz_h, *kk_horiz_d;
int ksize_horiz;

int  *bounds_vert_h, *bounds_vert_d;
double *kk_vert_h, *kk_vert_d;
int ksize_vert;

void coffes(unsigned int im_xsize, unsigned int im_ysize, int xsize, int ysize, unsigned char* data) {
    int *bounds_horiz;
    double *kk_horiz;

    float box[4] = {0, 0, 1.0f * im_xsize, 1.0f * im_ysize};

    precompute_coeffs(
            im_xsize,
            box[0],
            box[2],
            xsize,
            &bounds_horiz,
            &kk_horiz);

    double filterscale_horiz = (double)(box[2] - box[0]) / xsize;
    if (filterscale_horiz < 1.0) {
        filterscale_horiz = 1.0;
    }
    ksize_horiz = (int)ceil(filterscale_horiz) * 2 + 1;

//    normalize_coeffs_8bpc(xsize, ksize_horiz, kk_horiz);


    cudaMalloc(&kk_horiz_d, xsize * ksize_horiz * sizeof(double));
    cudaMalloc(&bounds_horiz_d, xsize * 2 * sizeof(int));

    coeffs<<<16, 256>>>(im_xsize,
                        box[0],
                        box[2],
                        xsize,
                        bounds_horiz_d,
                        kk_horiz_d);

    normalize_coeffs<<<xsize, ksize_horiz>>>(xsize, ksize_horiz, kk_horiz_d);

    bounds_horiz_h = (int*)malloc(xsize * 2 * sizeof(int));
    kk_horiz_h = (double*)malloc(xsize * ksize_horiz * sizeof(double));

    cudaMemcpy(bounds_horiz_h, bounds_horiz_d, xsize * 2 * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(kk_horiz_h, kk_horiz_d, xsize * ksize_horiz * sizeof(double), cudaMemcpyDeviceToHost);
//    for (int i = 0;i < )

    int *bounds_vert;
    double *kk_vert;

    precompute_coeffs(
            im_ysize,
            box[1],
            box[3],
            ysize,
            &bounds_vert,
            &kk_vert);

    double filterscale_vert = (double)(box[3] - box[1]) / ysize;
    if (filterscale_vert < 1.0) {
        filterscale_vert = 1.0;
    }
    ksize_vert = (int)ceil(filterscale_vert) * 2 + 1;
    cudaMalloc(&kk_vert_d, ysize * ksize_vert * sizeof(double));
    cudaMalloc(&bounds_vert_d, ysize * 2 * sizeof(int));

    normalize_coeffs_8bpc(ysize, ksize_vert, kk_vert);

    coeffs<<<16, 256>>>(im_ysize,
                        box[1],
                        box[3],
                        ysize,
                        bounds_vert_d,
                        kk_vert_d);

    normalize_coeffs<<<ysize, ksize_vert>>>(ysize, ksize_vert, kk_vert_d);

    bounds_vert_h = (int*)malloc(ysize * 2 * sizeof(int));
    kk_vert_h = (double*)malloc(ysize * ksize_vert * sizeof(double));

    cudaMemcpy(bounds_vert_h, bounds_vert_d, ysize * 2 * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(kk_vert_h, kk_vert_d, ysize * ksize_vert * sizeof(double), cudaMemcpyDeviceToHost);
    for (int i = 0;i < xsize * 2;i++) {
        if (bounds_horiz_h[i] != bounds_horiz[i]) {
            std::cout << "not equals! 11 " << std::endl;
        }
    }

//    for (int i = 0;i < xsize * ksize_horiz;i++) {
//        if (kk_horiz_h[i] != kk_horiz[i]) {
//            std::cout << "not equals! 12 " << std::endl;
//        }
//    }

    for (int i = 0;i < ysize * 2;i++) {
        if (bounds_vert_h[i] != bounds_vert[i]) {
            std::cout << "not equals! 21 " << std::endl;
        }
    }

    for (int i = 0;i < ysize * ksize_vert;i++) {
        if (kk_vert_h[i] != kk_vert[i]) {
            std::cout << "not equals! 22 " << std::endl;
        }
    }

//    // First used row in the source image
//    int ybox_first = bounds_vert[0];
//    // Last used row in the source image
    int ybox_last = bounds_vert[ysize * 2 - 2] + bounds_vert[ysize * 2 - 1];
    shift_ysize<<<16, 256>>>(bounds_vert_d, ysize);

    Imaging imagingInput = getImage(im_xsize, im_ysize, data);
    Imaging imagingTemp = getImage(xsize, im_ysize, NULL);
    ImagingResampleHorizontal_8bpc(imagingTemp, imagingInput, 0, ksize_horiz, bounds_horiz, kk_horiz);

    for (int i = 0;i < 20;i++) {
//        for (int j = 0;j < 3;j++) {
            UINT8 a1 = imagingTemp->image[0][i * 4 + 0];
            UINT8 a2 = imagingTemp->image[0][i * 4 + 1];
            UINT8 a3 = imagingTemp->image[0][i * 4 + 2];
            char output[1024];
            sprintf(output, "%d %d %d", (int)a1, (int)a2, (int)a3);
            std::cout<< output <<std::endl;
//        }
    }

}

int main() {
    const char* filename = "/home/siderzhang/2.jpg";
    UINT8 outBuffer[3600] = {
            144, 28, 63, 159, 41, 76, 159, 33, 71, 167, 36, 76, 175, 38, 82, 176, 36, 81, 184, 42, 88, 172, 30, 76, 184, 44, 89, 172, 32, 77, 182, 43, 85, 176, 35, 78, 187, 42, 85, 196, 46, 91, 177, 24, 70, 190, 40, 88, 162, 24, 75, 194, 56, 108, 173, 25, 77, 167, 15, 64, 180, 33, 75, 155, 21, 54, 142, 30, 52, 125, 29, 43, 119, 33, 42, 92, 9, 19, 108, 21, 37, 114, 29, 48, 92, 13, 34, 87, 19, 40, 33, 0, 3, 68, 17, 36, 88, 16, 36, 99, 15, 38, 117, 25, 48, 123, 28, 52, 112, 10, 34, 130, 26, 51, 104, 0, 23, 116, 12, 35, 110, 6, 39, 132, 27, 60, 141, 29, 65, 155, 39, 76, 158, 37, 78, 155, 30, 72, 159, 34, 76, 145, 23, 64, 136, 14, 55, 158, 38, 76, 177, 55, 94, 158, 36, 73, 149, 23, 60, 163, 33, 71, 169, 35, 72, 172, 40, 80, 167, 46, 91, 140, 19, 64, 163, 30, 77, 163, 26, 70, 147, 15, 54, 143, 24, 54, 124, 27, 46, 117, 35, 47, 91, 20, 26, 93, 23, 31, 94, 22, 33, 75, 4, 18, 67, 3, 19, 62, 10, 23, 24, 0, 0, 49, 11, 24, 82, 21, 37, 92, 19, 36, 109, 30, 49, 109, 24, 45, 108, 16, 37, 121, 28, 49, 100, 7, 26, 103, 10, 29, 67, 0, 14, 90, 10, 37, 102, 17, 46, 118, 31, 63, 124, 33, 68, 121, 30, 65, 133, 41, 78, 128, 38, 73, 129, 42, 76, 110, 25, 56, 100, 14, 43, 115, 29, 56, 122, 32, 57, 119, 25, 49, 126, 28, 51, 114, 19, 43, 116, 27, 57, 111, 22, 54, 122, 21, 55, 127, 22, 55, 115, 12, 41, 106, 17, 39, 81, 13, 24, 63, 11, 15, 70, 28, 29, 58, 20, 19, 82, 40, 42, 74, 32, 36, 63, 27, 31, 90, 64, 67, 106, 96, 95, 74, 60, 60, 69, 33, 37, 68, 17, 24, 79, 22, 31, 71, 8, 19, 84, 13, 27, 89, 17, 31, 84, 12, 26, 84, 12, 26, 44, 0, 3, 66, 1, 25, 77, 8, 36, 92, 23, 52, 88, 21, 52, 74, 9, 39, 76, 13, 44, 64, 5, 35, 64, 8, 35, 63, 9, 32, 56, 3, 23, 78, 24, 40, 75, 19, 32, 68, 8, 18, 95, 34, 42, 86, 26, 36, 104, 52, 65, 160, 105, 124, 81, 12, 33, 65, 0, 11, 79, 6, 26, 66, 4, 19, 64, 21, 28, 42, 14, 13, 21, 2, 0, 74, 56, 52, 100, 76, 74, 85, 61, 59, 113, 92, 89, 109, 96, 90, 116, 116, 106, 111, 108, 99, 121, 100, 97, 91, 59, 60, 80, 41, 44, 60, 18, 22, 73, 24, 30, 56, 3, 11, 60, 7, 15, 62, 9, 17, 70, 6, 32, 84, 21, 48, 87, 27, 55, 102, 46, 75, 107, 57, 86, 106, 61, 90, 121, 82, 111, 121, 87, 112, 107, 78, 100, 128, 101, 118, 86, 62, 75, 51, 26, 32, 63, 37, 40, 128, 100, 97, 188, 159, 155, 164, 139, 134, 124, 109, 106, 161, 142, 144, 70, 35, 42, 44, 1, 11, 45, 0, 10, 41, 4, 12, 30, 9, 14, 7, 0, 0, 50, 47, 42, 81, 76, 72, 149, 135, 132, 150, 131, 127, 109, 88, 83, 83, 69, 60, 119, 113, 101, 98, 96, 83, 99, 91, 80, 73, 58, 51, 76, 55, 50, 88, 63, 59, 106, 74, 75, 64, 30, 31, 47, 12, 16, 33, 0, 2, 89,
            12, 46, 124, 48, 84, 157, 87, 123, 197, 137, 171, 211, 162, 194, 205, 166, 197, 212, 181, 212, 206, 183, 209, 220, 205, 226, 221, 210, 226, 150, 143, 151, 85, 79, 81, 104, 96, 93, 163, 156, 146, 192, 184, 171, 171, 165, 151, 87, 84, 75, 37, 27, 25, 35, 8, 13, 50, 12, 23, 45, 5, 16, 113, 81, 92, 108, 91, 99, 111, 107, 108, 105, 107, 106, 134, 133, 129, 137, 126, 124, 121, 103, 99, 134, 114, 107, 110, 92, 82, 107, 95, 79, 154, 147, 129, 184, 178, 164, 148, 138, 128, 119, 102, 95, 126, 105, 102, 150, 122, 121, 141, 111, 113, 150, 119, 124, 145, 114, 120, 99, 0, 39, 152, 48, 97, 209, 116, 163, 253, 174, 219, 242, 180, 221, 193, 147, 184, 164, 130, 163, 139, 116, 145, 111, 97, 120, 147, 138, 155, 139, 134, 141, 94, 92, 93, 103, 102, 97, 127, 127, 115, 137, 138, 120, 174, 171, 156, 141, 130, 124, 54, 29, 32, 60, 17, 27, 73, 20, 36, 62, 7, 28, 162, 115, 135, 158, 130, 145, 172, 159, 168, 162, 157, 163, 141, 137, 138, 139, 125, 124, 133, 114, 110, 139, 114, 107, 150, 126, 114, 138, 121, 103, 179, 166, 147, 214, 198, 185, 208, 189, 182, 173, 148, 144, 182, 150, 151, 195, 155, 163, 201, 159, 169, 188, 146, 158, 150, 108, 122, 146, 7, 72, 170, 40, 104, 179, 64, 123, 176, 82, 135, 140, 66, 115, 100, 47, 89, 104, 68, 104, 108, 85, 114, 178, 164, 187, 197, 188, 205, 141, 138, 145, 36, 36, 36, 39, 40, 34, 84, 86, 73, 86, 89, 72, 128, 122, 108, 97, 69, 68, 85, 38, 48, 80, 12, 33, 112, 36, 62, 126, 50, 80, 192, 128, 155, 156, 114, 136, 135, 111, 125, 176, 166, 174, 201, 195, 197, 208, 197, 195, 191, 174, 167, 184, 162, 151, 169, 148, 131, 110, 95, 72, 185, 167, 147, 241, 211, 203, 237, 198, 199, 164, 117, 123, 145, 93, 105, 159, 101, 116, 220, 158, 179, 242, 180, 205, 213, 151, 176, 189, 26, 105, 220, 68, 143, 219, 85, 154, 174, 64, 125, 157, 74, 126, 157, 98, 142, 154, 117, 151, 214, 189, 218, 194, 179, 200, 184, 176, 189, 158, 151, 159, 90, 86, 87, 22, 18, 15, 49, 49, 41, 59, 59, 47, 105, 88, 81, 92, 41, 50, 113, 34, 56, 134, 39, 69, 168, 64, 99, 184, 84, 122, 236, 150, 185, 171, 112, 140, 167, 134, 151, 179, 166, 175, 198, 194, 195, 146, 141, 137, 139, 130, 121, 105, 92, 76, 111, 99, 77, 161, 154, 126, 193, 172, 151, 255, 218, 219, 200, 132, 145, 157, 80, 100, 152, 70, 93, 146, 57, 89, 156, 64, 101, 193, 101, 142, 221, 129, 170, 255, 83, 172, 233, 71, 154, 247, 109, 184, 215, 105, 168, 164, 84, 135, 190, 135, 174, 213, 180, 209, 206, 189, 208, 195, 183, 197, 166, 159, 167, 148, 138, 146, 120, 111, 114, 96, 86, 85, 24, 14, 12, 32, 25, 19, 82, 57, 60, 101, 44, 59, 142, 62
            , 89, 170, 70, 104, 183, 77, 115, 190, 85, 125, 247, 157, 193, 191, 131, 157, 180, 146, 162, 190, 179, 185, 185, 184, 182, 153, 150, 143, 167, 159, 148, 176, 160, 145, 179, 162, 142, 176, 162, 136, 216, 183, 168, 255, 183, 194, 198, 102, 129, 162, 58, 91, 175, 65, 104, 185, 68, 113, 188, 67, 120, 194, 72, 129, 195, 73, 132, 246, 77, 170, 250, 93, 180, 255, 125, 199, 209, 109, 169, 178, 113, 156, 198, 159, 188, 220, 201, 220, 222, 216, 226, 176, 171, 177, 129, 124, 128, 93, 82, 86, 129, 113, 116, 149, 130, 134, 48, 27, 32, 29, 8, 13, 46, 19, 28, 72, 34, 49, 114, 59, 82, 143, 68, 98, 162, 75, 109, 182, 95, 130, 242, 167, 198, 192, 143, 164, 172, 146, 157, 160, 154, 154, 152, 153, 147, 149, 145, 134, 160, 146, 135, 170, 144, 131, 171, 138, 123, 157, 122, 103, 216, 158, 154, 240, 142, 167, 197, 75, 116, 172, 42, 90, 189, 55, 108, 208, 66, 126, 211, 66, 133, 209, 63, 134, 204, 58, 131, 228, 69, 162, 253, 109, 196, 239, 123, 196, 188, 103, 158, 193, 143, 180, 198, 175, 195, 203, 196, 204, 211, 213, 212, 171, 170, 168, 137, 127, 126, 103, 83, 85, 158, 127, 135, 158, 120, 133, 131, 86, 106, 94, 49, 70, 80, 38, 60, 120, 89, 107, 123, 85, 108, 133, 73, 101, 153, 79, 112, 170, 93, 125, 206, 140, 167, 179, 138, 154, 187, 168, 172, 176, 175, 170, 175, 178, 167, 183, 177, 165, 179, 159, 148, 169, 132, 123, 182, 133, 126, 194, 142, 131, 216, 141, 146, 203, 89, 123, 192, 55, 107, 186, 42, 101, 204, 54, 118, 218, 61, 132, 221, 60, 138, 226, 64, 147, 235, 73, 158, 250, 106, 201, 234, 106, 191, 219, 119, 189, 187, 119, 168, 194, 160, 187, 191, 181, 190, 194, 196, 193, 157, 163, 153, 148, 144, 135, 130, 111, 107, 149, 110, 115, 162, 108, 122, 127, 59, 84, 151, 74, 108, 141, 58, 100, 120, 52, 91, 146, 112, 139, 100, 70, 96, 94, 39, 71, 130, 60, 94, 147, 74, 104, 158, 94, 118, 172, 134, 145, 236, 221, 218, 178, 180, 167, 165, 169, 154, 151, 144, 128, 151, 127, 115, 152, 107, 101, 197, 138, 134, 241, 173, 172, 173, 85, 99, 172, 51, 94, 187, 45, 105, 199, 50, 116, 216, 62, 134, 229, 68, 148, 225, 61, 147, 225, 61, 150, 233, 68, 160, 217, 86, 180, 215, 97, 183, 200, 114, 179, 170, 114, 159, 177, 153, 175, 142, 140, 143, 153, 161, 150, 148, 152, 137, 168, 156, 142, 141, 108, 103, 192, 132, 140, 173, 94, 116, 167, 69, 107, 186, 74, 124, 201, 83, 143, 179, 78, 134, 136, 75, 117, 61, 9, 47, 52, 0, 18, 122, 34, 76, 163, 73, 109, 147, 71, 97, 155, 108, 118, 214, 193, 188, 174, 175, 159, 151, 159, 138, 120, 117, 98, 110, 93, 77, 121, 84, 75, 185, 134, 130, 255, 198, 197, 178, 95, 111, 183, 66, 111, 200, 59, 125, 200, 53, 123, 209, 57, 132, 229, 69, 155, 233, 70, 161, 227, 63, 158, 232, 68, 165, 177, 52, 146, 188, 78, 164, 196, 115, 181, 169, 118, 161, 149, 130, 149, 89, 88, 86, 141, 145, 130, 181, 178, 159, 172, 149, 135, 168, 117, 114, 179, 97, 111, 170, 61, 92, 205, 73, 123, 215, 69, 134, 229, 72, 149, 213, 71, 147, 195, 85, 148, 135, 31, 90, 130, 8, 69, 190, 58, 116, 211, 83, 132, 155, 50, 83, 129, 61, 76, 154, 119, 117, 154, 151, 132, 138, 151, 125, 121, 131, 106, 91, 93, 71, 100, 87, 70, 133, 110, 96, 191, 161, 151, 176, 119, 128, 185, 79, 125, 208, 75, 141, 200, 60, 131, 201, 54, 132, 224, 70, 158, 234, 76, 169, 230, 72, 169, 237, 78, 178, 207, 87, 176, 154, 51, 130, 199, 134, 190, 200, 163, 197, 124, 108, 119, 80, 77, 72, 190, 188, 173, 191, 182, 165, 176, 148, 137, 227, 173, 173, 183, 94, 114, 189, 70, 108, 214, 70, 129, 234, 72, 145, 241, 65, 150, 250, 82, 167, 231, 90, 168, 215, 81, 154, 224, 77, 148, 240, 91, 157, 211, 69, 125, 134, 18, 57, 125, 46, 65, 158, 116, 117, 157, 150, 132, 150, 160, 135, 158, 173, 144, 130, 140, 113, 143, 142, 121, 129, 118, 100, 131, 112, 98, 157, 115, 119, 140, 58, 94, 182, 75, 129, 195, 75, 137, 202, 71, 141, 226, 80, 163, 231, 78, 168, 227, 69, 166, 234, 75, 175, 237, 120, 198, 200, 108, 171, 189, 147, 183, 183, 169, 186, 132, 120, 124, 130, 119, 113, 186, 173, 164, 197, 180, 170, 229, 205, 201, 213, 168, 175, 154, 74, 101, 184, 74, 119, 209, 77, 137, 231, 77, 149, 251, 73, 157, 255, 80, 168, 230, 74, 158, 251, 105, 186, 249, 101, 177, 203, 59, 128, 170, 38, 95, 173, 64, 105, 173, 95, 117, 153, 108, 113, 126, 112, 101, 105, 110, 88, 115, 129, 103, 151, 164, 138, 144, 142, 121, 112, 96, 80, 114, 84, 74, 116, 76, 76, 157, 109, 123, 198, 135, 162, 155, 70, 109, 187, 80, 134, 211, 77, 148, 230, 78, 161, 249, 85, 180, 230, 61, 162, 223, 116, 186, 203, 124, 179, 181, 161, 186, 141, 144, 151, 95, 91, 90, 141, 130, 124, 216, 203, 195, 210, 195, 188, 217, 199, 199, 141, 105, 117, 131, 62, 93, 165, 67, 114, 196, 78, 138, 222, 80, 152, 250, 81, 162, 251, 77, 164, 255, 102, 190, 240, 88, 173, 217, 69, 147, 198, 57, 126, 185, 56, 113, 181, 74, 118, 193, 113, 140, 201, 150, 159, 175, 154, 149, 170, 167, 152, 155, 160, 140, 119, 122, 101, 127, 118, 103, 136, 112, 102, 119, 78, 74, 135, 91, 92, 181, 146, 152, 212, 170, 184, 182, 116, 144, 175, 79, 125, 207, 81, 145, 221, 72, 153, 223, 60, 153, 238, 68, 167, 164, 75, 141, 203, 143, 194, 167, 162, 185, 104, 122, 126, 108, 114, 112, 142, 142, 134, 172, 165, 157, 195
            , 186, 181, 136, 126, 127, 86, 60, 73, 142, 83, 113, 124, 41, 87, 115, 12, 68, 165, 40, 108, 231, 83, 159, 234, 75, 159, 237, 83, 169, 213, 61, 146, 200, 51, 132, 211, 67, 139, 215, 83, 143, 203, 91, 139, 191, 106, 137, 185, 127, 142, 180, 148, 151, 165, 148, 141, 154, 145, 136, 125, 112, 103, 145, 120, 116, 131, 91, 92, 98, 43, 49, 182, 126, 137, 192, 152, 163, 157, 112, 132, 180, 110, 144, 184, 88, 136, 185, 58, 125, 213, 66, 147, 234, 75, 168, 240, 74, 174, 161, 82, 147, 160, 108, 157, 155, 157, 180, 128, 151, 157, 108, 116, 118, 123, 119, 118, 141, 131, 129, 147, 135, 135, 72, 59, 66, 51, 22, 40, 95, 38, 73, 83, 4, 51, 75, 0, 38, 139, 26, 92, 220, 84, 160, 232, 85, 166, 201, 52, 136, 215, 63, 148, 224, 75, 156, 228, 83, 158, 233, 98, 164, 224, 107, 159, 187, 93, 130, 141, 71, 97, 127, 79, 93, 118, 85, 92, 141, 112, 116, 157, 126, 131, 167, 124, 133, 118, 60, 75, 95, 22, 42, 194, 122, 144, 216, 164, 186, 142, 89, 115, 133, 58, 97, 162, 63, 117, 175, 48, 119, 224, 78, 161, 241, 83, 178, 236, 74, 173, 178, 103, 168, 181, 133, 183, 138, 140, 165, 98, 116, 128, 122, 122, 130, 114, 101, 108, 98, 75, 83, 158, 130, 142, 217, 188, 206, 140, 100, 127, 79, 12, 55, 130, 45, 100, 176, 82, 142, 201, 92, 159, 205, 77, 152, 194, 52, 134, 219, 71, 157, 245, 95, 182, 255, 104, 187, 239, 92, 170, 231, 94, 162, 230, 106, 166, 200, 94, 142, 154, 68, 105, 108, 42, 70, 131, 76, 99, 155, 103, 125, 166, 111, 134, 161, 92, 120, 143, 62, 95, 150, 54, 92, 162, 70, 109, 238, 173, 205, 241, 180, 214, 154, 72, 118, 145, 41, 100, 192, 64, 138, 234, 90, 175, 214, 59, 153, 230, 70, 170, 211, 130, 199, 163, 105, 163, 150, 140, 174, 119, 121, 142, 97, 78, 98, 123, 88, 108, 138, 93, 114, 149, 96, 122, 142, 88, 120, 165, 103, 144, 139, 54, 109, 182, 83, 147, 186, 81, 148, 185, 69, 142, 212, 79, 158, 255, 112, 197, 243, 95, 183, 254, 104, 193, 255, 106, 191, 245, 98, 178, 233, 91, 165, 223, 92, 158, 206, 91, 148, 186, 86, 136, 160, 75, 117, 163, 87, 125, 165, 91, 128, 168, 87, 128, 160, 69, 113, 179, 75, 126, 208, 89, 145, 169, 58, 111, 200, 119, 162, 255, 207, 250, 251, 158, 211, 205, 95, 158, 197, 66, 142, 221, 74, 162, 219, 65, 161, 227, 68, 168, 215, 114, 192, 179, 103, 168, 175, 144, 188, 146, 125, 158, 136, 89, 123, 156, 90, 126, 156, 75, 116, 173, 86, 131, 163, 77, 126, 194, 99, 157, 194, 80, 150, 199, 74, 152, 218, 91, 172, 234, 101, 184, 239, 92, 180, 241, 91, 180, 249, 99, 188, 248, 98, 187, 249, 97, 184, 246, 97, 181, 240, 95, 174, 229, 92, 164, 217, 90, 157, 207, 91, 152, 213, 108, 164, 184, 85, 139, 189, 90, 145, 196, 91, 148, 195, 79, 142, 214, 87, 156, 227, 90, 162, 217, 87, 157, 183, 84, 139, 214, 124, 178, 255, 170, 232, 235, 114, 185, 208, 71, 153, 216, 67, 157, 243, 87, 186, 233, 74, 176, 241, 127, 212, 174, 82, 157, 200, 152, 204, 198, 157, 201, 159, 90, 137, 178, 86, 137, 194, 87, 143, 199, 85, 145, 209, 96, 160, 199, 79, 151, 239, 101, 186, 228, 84, 172, 246, 106, 195, 250, 106, 194, 252, 97, 189, 255, 107, 198, 255, 108, 199, 255, 107, 197, 250, 97, 187, 240, 90, 177, 243, 96, 177, 250, 107, 186, 243, 107, 183, 226, 99, 170, 225, 104, 173, 210, 91, 159, 241, 122, 190, 230, 104, 177, 229, 94, 171, 247, 104, 186, 232, 78, 164, 246, 100, 181, 238, 121, 189, 170, 61, 126, 213, 93, 164, 177, 44, 123, 244, 100, 187, 237, 83, 177, 233, 74, 174, 242, 83, 185, 242, 117, 208, 210, 107, 186, 219, 158, 218, 181, 124, 177, 159, 71, 129, 170, 60, 121, 216, 89, 156, 218, 82, 156, 234, 102, 178, 241, 101, 187, 254, 100, 194, 255, 103, 200, 255, 105, 200, 254, 103, 198, 255, 97, 190, 254, 95, 188, 246, 91, 183, 245, 95, 185, 250, 97, 187, 252, 100, 187, 250, 100, 187, 247, 99, 183, 243, 100, 182, 238, 99, 180, 235, 100, 179, 236, 101, 180, 236, 99, 181, 236, 94, 180, 241, 93, 181, 251, 93, 186, 255, 95, 190, 255, 100, 191, 218, 86, 162, 199, 73, 147, 203, 70, 149, 212, 70, 156, 220, 69, 162, 253, 95, 192, 217, 56, 159, 237, 76, 180, 245, 118, 209, 213, 106, 188, 226, 160, 222, 209, 146, 203, 177, 83, 145, 189, 70, 138, 232, 96, 172, 242, 99, 179, 250, 108, 194, 254, 107, 198, 255, 104, 204, 255, 103, 206, 255, 105, 202, 251, 102, 196, 253, 95, 190, 251, 93, 186, 251, 98, 189, 252, 102, 192, 255, 102, 192, 252, 99, 189, 245, 95, 184, 242, 92, 181, 242, 94, 182, 245, 97, 185, 244, 97, 185, 241, 94, 182, 241, 92, 182, 245, 91, 185, 252, 94, 191, 255, 93, 192, 255, 87, 188, 245, 84, 180, 215, 69, 154, 223, 84, 165, 223, 79, 166, 211, 62, 152, 210, 54, 151, 255, 96, 196, 237, 76, 180, 221, 60, 164, 236, 113, 203, 200, 97, 178, 208, 146, 209, 207, 148, 206, 164, 72, 137, 186, 69, 140, 229, 94, 173, 249, 107, 193, 249, 110, 201, 250, 109, 205, 255, 104, 205, 255, 103, 204, 248, 104, 201, 243, 102, 196, 245, 98, 189, 246, 96, 186, 250, 100, 189, 255, 102, 192, 255, 103, 193, 250, 97, 187, 243, 90, 180, 241, 88, 179, 246, 93, 184, 252, 99, 190, 248, 94, 188, 244, 89, 183, 244, 86, 183, 251, 91, 189, 255, 97, 198, 255, 93, 195, 247, 81, 183, 234, 70, 169, 226, 73, 163, 242, 92, 179, 230, 80, 169, 217, 64, 155, 200, 44, 141, 238, 79, 179, 254, 95
            , 198, 222, 62, 168, 228, 105, 195, 194, 91, 172, 184, 122, 185, 181, 122, 180, 140, 49, 116, 184, 68, 143, 222, 89, 172, 238, 98, 187, 241, 104, 196, 241, 101, 198, 246, 97, 201, 245, 96, 198, 238, 100, 196, 236, 101, 193, 241, 101, 190, 246, 99, 188, 243, 93, 182, 246, 93, 183, 246, 93, 183, 245, 92, 182, 244, 91, 182, 245, 92, 183, 250, 95, 189, 253, 98, 192, 249, 91, 188, 245, 85, 183, 246, 84, 185, 253, 91, 192, 255, 95, 198, 251, 89, 192, 239, 77, 180, 229, 69, 167, 251, 93, 186, 243, 88, 179, 220, 65, 157, 229, 75, 169, 202, 48, 144, 202, 47, 148, 251, 96, 198, 237, 82, 186, 226, 99, 190, 204, 97, 179, 181, 115, 179, 166, 105, 164, 149, 57, 124, 211, 93, 169, 234, 100, 185, 235, 94, 186, 234, 99, 194, 235, 96, 197, 239, 92, 196, 239, 92, 196, 235, 96, 195, 235, 100, 194, 241, 101, 190, 248, 101, 190, 240, 90, 179, 238, 86, 173, 236, 83, 173, 238, 88, 177, 245, 95, 185, 251, 98, 189, 249, 94, 188, 247, 89, 184, 245, 85, 183, 245, 83, 182, 248, 86, 187, 251, 91, 191, 249, 89, 189, 238, 82, 181, 234, 78, 177, 235, 79, 176, 252, 94, 187, 236, 79, 170, 212, 59, 150, 233, 82, 175, 204, 56, 150, 174, 28, 126, 233, 87, 186, 231, 85, 186, 220, 91, 183, 212, 105, 187, 187, 118, 183, 162, 101, 160, 185, 89, 160, 242, 124, 202, 248, 114, 199, 239, 98, 190, 238, 102, 200, 239, 100, 201, 242, 95, 201, 240, 93, 197, 233, 97, 195, 234, 101, 194, 241, 101, 190, 245, 101, 188, 245, 92, 182, 238, 86, 173, 232, 79, 169, 235, 85, 174, 244, 94, 184, 248, 98, 188, 245, 90, 184, 238, 80, 175, 243, 83, 181, 245, 83, 182, 249, 87, 188, 251, 91, 191, 243, 84, 184, 230, 76, 174, 230, 79, 176, 240, 89, 184, 231, 73, 166, 232, 75, 166, 219, 66, 157, 225, 76, 168, 204, 59, 152, 172, 28, 125, 217, 75, 173, 211, 69, 169};
    unsigned int im_xsize = 40;
    unsigned int im_ysize = 30;

//    unsigned int width;
//    unsigned int height;
    unsigned int channels = 3;

//    read_jpeg_file(filename, &outBuffer, &im_xsize, &im_ysize, &channels);
//    free(outBuffer)


    for (int j = 0;j < 3;j++) {
        for (int i = 0;i < 1200;i++) {
            uint s = ((uint8_t *)outBuffer)[i * 3 + j];
            char buffer[10];
            sprintf(buffer, "%d", s);
            std::cout<<buffer<<", ";
        }
        std::cout<<std::endl;
    }

    int xsize = 20;

    int ysize = 15;
    coffes(im_xsize, im_ysize, xsize, ysize, outBuffer);

//    cudaMalloc(&kk_vert_d, ysize * ksize_vert * sizeof(double));
//    cudaMalloc(&bounds_vert_d, ysize * 2 * sizeof(int));
    unsigned char* input;//  (unsigned char*)calloc(im_xsize * im_ysize, 3 * sizeof(unsigned char));
    cudaMalloc(&input, im_xsize * im_ysize * 3 * sizeof(unsigned char));

    cudaMemcpy(input, outBuffer, xsize * im_ysize * 3 * sizeof(unsigned char), cudaMemcpyHostToDevice);

    unsigned char* temp;// = (unsigned char*)calloc(xsize * ysize, 3 * sizeof(unsigned char));
    cudaMalloc(&temp, xsize * im_ysize * 3 * sizeof(unsigned char));

    unsigned char* output;// = (unsigned char*)calloc(xsize * ysize, 3 * sizeof(unsigned char));
    cudaMalloc(&output, xsize * ysize * 3 * sizeof(unsigned char));

    unsigned char* looksup_d;
    cudaMalloc(&looksup_d, 1280 * sizeof(unsigned char));
    cudaMemcpy(looksup_d, lookups_h, 1280 * sizeof(unsigned char), cudaMemcpyHostToDevice);

    build_result_horiz<<<16, 128>>>(im_xsize, 3, xsize, im_ysize, input, temp, ksize_horiz, bounds_horiz_d, kk_horiz_d, looksup_d);
    build_result_vert<<<16, 128>>>(xsize, 3, xsize, ysize, temp, output, ksize_vert, bounds_vert_d, kk_vert_d, looksup_d);

    unsigned char* tempTest = (unsigned char*)calloc(xsize * ysize, 3 * sizeof(unsigned char));
    cudaMemcpy(tempTest, output, xsize * ysize * 3 * sizeof(unsigned char), cudaMemcpyDeviceToHost);

    for (int j = 0;j < 3;j++) {
        for (int i = 0;i < 300;i++) {
            uint s = ((uint8_t *)tempTest)[i * 3 + j];
            char buffer[10];
            sprintf(buffer, "%d", s);
            std::cout<<buffer<<", ";
        }
        std::cout<<std::endl;
    }
//    build_result_vert<<<16, 256>>>(im_xsize, 3, xsize, ysize, temp, output, ksize_vert, bounds_vert_d, kk_vert_d);

    cudaDeviceSynchronize();

    std::cout<<"Hello World!"<<std::endl;
    return 0;
}