//
// Created by siderzhangPC on 2024/4/19.
//

#ifndef CUDAHELLOWORLD_IPL_H
#define CUDAHELLOWORLD_IPL_H

#include "Imaging.h"

struct filter {
    double (*filter)(double x);
    double support;
};
int
precompute_coeffs(
        int inSize,
        float in0,
        float in1,
        int outSize,
//        struct filter *filterp,
        int **boundsp,
        double **kkp);
void
normalize_coeffs_8bpc(int outSize, int ksize, double *prekk);

Imaging getImage(int width, int height, unsigned char* data);

void
ImagingResampleHorizontal_8bpc(
        Imaging imOut, Imaging imIn, int offset, int ksize, int *bounds, double *prekk);

#endif //CUDAHELLOWORLD_IPL_H
