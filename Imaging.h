//
// Created by siderzhangPC on 2024/4/20.
//

#ifndef CUDAHELLOWORLD_IMAGING_H
#define CUDAHELLOWORLD_IMAGING_H

typedef struct ImagingMemoryInstance *Imaging;
typedef struct ImagingPaletteInstance *ImagingPalette;
#define IMAGING_TYPE_UINT8 0
#define IMAGING_TYPE_INT32 1
#define IMAGING_TYPE_FLOAT32 2
#define IMAGING_TYPE_SPECIAL 3 /* check mode for details */

typedef signed char         INT8, *PINT8;
typedef signed short        INT16, *PINT16;
typedef signed int          INT32, *PINT32;
typedef unsigned char       UINT8, *PUINT8;
typedef unsigned short      UINT16, *PUINT16;
typedef unsigned int        UINT32, *PUINT32;

#define PRECISION_BITS (32 - 8 - 2)

#define IMAGING_MODE_LENGTH \
    6 + 1 /* Band names ("1", "L", "P", "RGB", "RGBA", "CMYK", "YCbCr", "BGR;xy") */

typedef struct {
    char *ptr;
    int size;
} ImagingMemoryBlock;


struct ImagingMemoryInstance {
    /* Format */
    char mode[IMAGING_MODE_LENGTH]; /* Band names ("1", "L", "P", "RGB", "RGBA", "CMYK",
                                       "YCbCr", "BGR;xy") */
    int type;                       /* Data type (IMAGING_TYPE_*) */
    int depth;                      /* Depth (ignored in this version) */
    int bands;                      /* Number of bands (1, 2, 3, or 4) */
    int xsize;                      /* Image dimension. */
    int ysize;

    /* Colour palette (for "P" images only) */
    ImagingPalette palette;

    /* Data pointers */
    UINT8 **image8;  /* Set for 8-bit images (pixelsize=1). */
    INT32 **image32; /* Set for 32-bit images (pixelsize=4). */

    /* Internals */
    char **image;               /* Actual raster data. */
    char *block;                /* Set if data is allocated in a single block. */
    ImagingMemoryBlock *blocks; /* Memory blocks for pixel storage */

    int pixelsize; /* Size of a pixel, in bytes (1, 2 or 4) */
    int linesize;  /* Size of a line, in bytes (xsize * pixelsize) */

    /* Virtual methods */
    void (*destroy)(Imaging im);
};


struct ImagingPaletteInstance {
    /* Format */
    char mode[IMAGING_MODE_LENGTH]; /* Band names */

    /* Data */
    int size;
    UINT8 palette[1024]; /* Palette data (same format as image data) */

    INT16 *cache;   /* Palette cache (used for predefined palettes) */
    int keep_cache; /* This palette will be reused; keep cache */
};

#endif //CUDAHELLOWORLD_IMAGING_H
