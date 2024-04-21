//
// Created by siderzhangPC on 2024/4/19.
//

#include "ipl.h"
#include "Imaging.h"
#include <math.h>
#include <cstdio>


#define MAKE_UINT32(u0, u1, u2, u3) \
    ((UINT32)(u0) | ((UINT32)(u1) << 8) | ((UINT32)(u2) << 16) | ((UINT32)(u3) << 24))

/* Handles values form -640 to 639. */
UINT8 _clip8_lookups[1280] = {
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

UINT8 *clip8_lookups = &_clip8_lookups[640];

static inline UINT8
clip8(int in) {
    return clip8_lookups[in >> PRECISION_BITS];
}

static inline double
bilinear_filter(double x) {
    if (x < 0.0) {
        x = -x;
    }
    if (x < 1.0) {
        return 1.0 - x;
    }
    return 0.0;
}

static struct filter BILINEAR = {bilinear_filter, 1.0};

int
precompute_coeffs(
        int inSize,
        float in0,
        float in1,
        int outSize,
//        struct filter *filterp,
        int **boundsp,
        double **kkp) {
    double support, scale, filterscale;
    double center, ww, ss;
    int xx, x, ksize, xmin, xmax;
    int *bounds;
    double *kk, *k;

    // panding
    struct filter *filterp = &BILINEAR;

    /* prepare for horizontal stretch */
    filterscale = scale = (double)(in1 - in0) / outSize;
    if (filterscale < 1.0) {
        filterscale = 1.0;
    }

    /* determine support size (length of resampling filter) */
    support = filterp->support * filterscale;

    /* maximum number of coeffs */
    ksize = (int)ceil(support) * 2 + 1;

    // check for overflow
    if (outSize > INT_MAX / (ksize * (int)sizeof(double))) {
        return 0;
    }

    /* coefficient buffer */
    /* malloc check ok, overflow checked above */
    kk = (double*)malloc(outSize * ksize * sizeof(double));
    if (!kk) {
        return 0;
    }

    /* malloc check ok, ksize*sizeof(double) > 2*sizeof(int) */
    bounds = (int*)malloc(outSize * 2 * sizeof(int));
    if (!bounds) {
        free(kk);
        return 0;
    }

    for (xx = 0; xx < outSize; xx++) {
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
            double w = filterp->filter((x + xmin - center + 0.5) * ss);
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
    *boundsp = bounds;
    *kkp = kk;
    return ksize;
}


void
normalize_coeffs_8bpc(int outSize, int ksize, double *prekk) {
    int x;
    INT32 *kk;

    // use the same buffer for normalized coefficients
    kk = (INT32 *)prekk;

    for (x = 0; x < outSize * ksize; x++) {

        if (x == 15) {
            int sd = 0;
        }
        if (prekk[x] < 0) {
            kk[x] = (int)(-0.5 + prekk[x] * (1 << PRECISION_BITS));
        } else {
            kk[x] = (int)(0.5 + prekk[x] * (1 << PRECISION_BITS));
        }
    }
}

void
ImagingResampleHorizontal_8bpc(
        Imaging imOut, Imaging imIn, int offset, int ksize, int *bounds, double *prekk) {
//    ImagingSectionCookie cookie;
    int ss0, ss1, ss2, ss3;
    int xx, yy, x, xmin, xmax;
    INT32 *k, *kk;

    // use the same buffer for normalized coefficients
    kk = (INT32 *)prekk;
    normalize_coeffs_8bpc(imOut->xsize, ksize, prekk);

//    ImagingSectionEnter(&cookie);
    if (imIn->image8) {
        for (yy = 0; yy < imOut->ysize; yy++) {
            for (xx = 0; xx < imOut->xsize; xx++) {
                xmin = bounds[xx * 2 + 0];
                xmax = bounds[xx * 2 + 1];
                k = &kk[xx * ksize];
                ss0 = 1 << (PRECISION_BITS - 1);
                for (x = 0; x < xmax; x++) {
                    ss0 += ((UINT8)imIn->image8[yy + offset][x + xmin]) * k[x];
                }
                imOut->image8[yy][xx] = clip8(ss0);
            }
        }
    } else if (imIn->type == IMAGING_TYPE_UINT8) {
        if (imIn->bands == 2) {
            for (yy = 0; yy < imOut->ysize; yy++) {
                for (xx = 0; xx < imOut->xsize; xx++) {
                    UINT32 v;
                    xmin = bounds[xx * 2 + 0];
                    xmax = bounds[xx * 2 + 1];
                    k = &kk[xx * ksize];
                    ss0 = ss3 = 1 << (PRECISION_BITS - 1);
                    for (x = 0; x < xmax; x++) {
                        ss0 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 0]) *
                               k[x];
                        ss3 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 3]) *
                               k[x];
                    }
                    v = MAKE_UINT32(clip8(ss0), 0, 0, clip8(ss3));
                    memcpy(imOut->image[yy] + xx * sizeof(v), &v, sizeof(v));
                }
            }
        } else if (imIn->bands == 3) {
            for (yy = 0; yy < imOut->ysize; yy++) {
                for (xx = 0; xx < imOut->xsize; xx++) {
                    UINT32 v;
                    xmin = bounds[xx * 2 + 0];
                    xmax = bounds[xx * 2 + 1];
                    k = &kk[xx * ksize];
                    ss0 = ss1 = ss2 = 1 << (PRECISION_BITS - 1);
                    for (x = 0; x < xmax; x++) {
                        ss0 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 0]) *
                               k[x];
                        ss1 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 1]) *
                               k[x];
                        ss2 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 2]) *
                               k[x];
                    }
                    UINT8 as = clip8(ss0);
                    v = MAKE_UINT32(clip8(ss0), clip8(ss1), clip8(ss2), 0);
                    memcpy(imOut->image[yy] + xx * sizeof(v), &v, sizeof(v));
                }
            }
        } else {
            for (yy = 0; yy < imOut->ysize; yy++) {
                for (xx = 0; xx < imOut->xsize; xx++) {
                    UINT32 v;
                    xmin = bounds[xx * 2 + 0];
                    xmax = bounds[xx * 2 + 1];
                    k = &kk[xx * ksize];
                    ss0 = ss1 = ss2 = ss3 = 1 << (PRECISION_BITS - 1);
                    for (x = 0; x < xmax; x++) {
                        ss0 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 0]) *
                               k[x];
                        ss1 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 1]) *
                               k[x];
                        ss2 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 2]) *
                               k[x];
                        ss3 += ((UINT8)imIn->image[yy + offset][(x + xmin) * 4 + 3]) *
                               k[x];
                    }
                    v = MAKE_UINT32(clip8(ss0), clip8(ss1), clip8(ss2), clip8(ss3));
                    memcpy(imOut->image[yy] + xx * sizeof(v), &v, sizeof(v));
                }
            }
        }
    }
//    ImagingSectionLeave(&cookie);
}

void
ImagingResampleVertical_8bpc(
        Imaging imOut, Imaging imIn, int offset, int ksize, int *bounds, double *prekk) {
//    ImagingSectionCookie cookie;
    int ss0, ss1, ss2, ss3;
    int xx, yy, y, ymin, ymax;
    INT32 *k, *kk;

    // use the same buffer for normalized coefficients
    kk = (INT32 *)prekk;
    normalize_coeffs_8bpc(imOut->ysize, ksize, prekk);

//    ImagingSectionEnter(&cookie);
    if (imIn->image8) {
        for (yy = 0; yy < imOut->ysize; yy++) {
            k = &kk[yy * ksize];
            ymin = bounds[yy * 2 + 0];
            ymax = bounds[yy * 2 + 1];
            for (xx = 0; xx < imOut->xsize; xx++) {
                ss0 = 1 << (PRECISION_BITS - 1);
                for (y = 0; y < ymax; y++) {
                    ss0 += ((UINT8)imIn->image8[y + ymin][xx]) * k[y];
                }
                imOut->image8[yy][xx] = clip8(ss0);
            }
        }
    } else if (imIn->type == IMAGING_TYPE_UINT8) {
        if (imIn->bands == 2) {
            for (yy = 0; yy < imOut->ysize; yy++) {
                k = &kk[yy * ksize];
                ymin = bounds[yy * 2 + 0];
                ymax = bounds[yy * 2 + 1];
                for (xx = 0; xx < imOut->xsize; xx++) {
                    UINT32 v;
                    ss0 = ss3 = 1 << (PRECISION_BITS - 1);
                    for (y = 0; y < ymax; y++) {
                        ss0 += ((UINT8)imIn->image[y + ymin][xx * 4 + 0]) * k[y];
                        ss3 += ((UINT8)imIn->image[y + ymin][xx * 4 + 3]) * k[y];
                    }
                    v = MAKE_UINT32(clip8(ss0), 0, 0, clip8(ss3));
                    memcpy(imOut->image[yy] + xx * sizeof(v), &v, sizeof(v));
                }
            }
        } else if (imIn->bands == 3) {
            for (yy = 0; yy < imOut->ysize; yy++) {
                k = &kk[yy * ksize];
                ymin = bounds[yy * 2 + 0];
                ymax = bounds[yy * 2 + 1];
                for (xx = 0; xx < imOut->xsize; xx++) {
                    UINT32 v;
                    ss0 = ss1 = ss2 = 1 << (PRECISION_BITS - 1);
                    for (y = 0; y < ymax; y++) {
                        ss0 += ((UINT8)imIn->image[y + ymin][xx * 4 + 0]) * k[y];
                        ss1 += ((UINT8)imIn->image[y + ymin][xx * 4 + 1]) * k[y];
                        ss2 += ((UINT8)imIn->image[y + ymin][xx * 4 + 2]) * k[y];
                    }
                    v = MAKE_UINT32(clip8(ss0), clip8(ss1), clip8(ss2), 0);
                    memcpy(imOut->image[yy] + xx * sizeof(v), &v, sizeof(v));
                }
            }
        } else {
            for (yy = 0; yy < imOut->ysize; yy++) {
                k = &kk[yy * ksize];
                ymin = bounds[yy * 2 + 0];
                ymax = bounds[yy * 2 + 1];
                for (xx = 0; xx < imOut->xsize; xx++) {
                    UINT32 v;
                    ss0 = ss1 = ss2 = ss3 = 1 << (PRECISION_BITS - 1);
                    for (y = 0; y < ymax; y++) {
                        ss0 += ((UINT8)imIn->image[y + ymin][xx * 4 + 0]) * k[y];
                        ss1 += ((UINT8)imIn->image[y + ymin][xx * 4 + 1]) * k[y];
                        ss2 += ((UINT8)imIn->image[y + ymin][xx * 4 + 2]) * k[y];
                        ss3 += ((UINT8)imIn->image[y + ymin][xx * 4 + 3]) * k[y];
                    }
                    v = MAKE_UINT32(clip8(ss0), clip8(ss1), clip8(ss2), clip8(ss3));
                    memcpy(imOut->image[yy] + xx * sizeof(v), &v, sizeof(v));
                }
            }
        }
    }
//    ImagingSectionLeave(&cookie);
}

Imaging getImage(int width, int height, unsigned char* data) {
    Imaging imaging = (Imaging)calloc(1, sizeof(ImagingMemoryInstance));
    memcpy(imaging->mode, "RGB", sizeof("RGB"));
    imaging->bands = 3;
    imaging->xsize = width;
    imaging->ysize = height;
    imaging->image = static_cast<char **>(calloc(imaging->ysize, sizeof(char *)));
    for (int i = 0;i < imaging->ysize;i++) {
        imaging->image[i] = (char*)calloc(imaging->xsize * 4, sizeof(UINT8));
        if (data != NULL) {
            for (int j = 0;j < imaging->xsize;j++) {
                imaging->image[i][j * 4 + 0] = data[i * imaging->xsize * 3 + j * 3 + 0];
                imaging->image[i][j * 4 + 1] = data[i * imaging->xsize * 3 + j * 3 + 1];
                imaging->image[i][j * 4 + 2] = data[i * imaging->xsize * 3 + j * 3 + 2];
                imaging->image[i][j * 4 + 3] = 0;
            }
        }
    }

    return imaging;
}