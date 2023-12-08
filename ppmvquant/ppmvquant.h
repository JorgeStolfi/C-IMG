/* ppmvquant.h - options, types and prototypes for ppmvquant.c */
/* Last edited on 2017-06-20 20:43:49 by stolfilocal */

#ifndef ppmvquant_H
#define ppmvquant_H

#include <uint16_image_RGB_hist.h>
#include <uint16_image_RGB_table.h>
#include <uint16_image_RGB_medcut.h>
#include <yuvhacks.h>
#include <uint16_image.h>

typedef struct rgb_pixel_t { uint16_t c[3]; } rgb_pixel_t;

#define MAXCOLORS 32767

/* Options for the median cut algorithm: */

/* How to choose the axis of greatest variation: */
#define LARGE_LUM
/* #define LARGE_NORM */

/* How to choose the representative color in each box: */
#define REP_AVERAGE_PIXELS
/* #define REP_CENTER_BOX */
/* #define REP_AVERAGE_COLORS */

typedef int matchfn_t
  ( long r, long g, long b,
    int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors, 
    int mapmaxval
  );
  /* A function that finds the best match to (r,g,b)/maxval in the 
    colormap cm, whose colors range in [0..mapmaxval].
    Note that (r,g,b) may be somewhat outside of the color cube. */

typedef struct box
  { int ind;
    int colors;
    long sum;
  } box;
  
typedef box* box_vector;

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { bool_t floyd;
    matchfn_t *match;
    char *mapname; 
    int ncolors; 
    char *imgname;
    bool_t verbose;
  }options_t;

/* PROTOTYPES */

int main(int argc, char**argv);

options_t *parse_options(int argc, char **argv);

uint16_image_RGB_hist_vector choose_colormap(uint16_image_t *img, int *newcolorsp);

uint16_image_RGB_hist_vector collect_colors(uint16_image_t *map, int *newcolorsp);

void uint16_image_quantize_floyd
  ( uint16_image_t *img, 
    uint16_image_RGB_hist_vector cm, 
    int newcolors,
    int mapmaxval,
    matchfn_t *match
  );
  /* Replaces each pixel in the image by a pixel from the given colormap. */

void simple_quantize_uint16_image
  ( uint16_image_t *img, 
    uint16_image_RGB_hist_vector cm, 
    int newcolors,
    int mapmaxval,
    matchfn_t *match
  );

ppm_pixel_t *choose_color
  ( long r, long g, long b, int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval,
    matchfn_t *match,
    uint16_image_RGB_table cht, 
    int *addtohash
  );
  /* Finds best match to (r,g,b) in colormap cm, using/updating the 
    hash table if appropriate. */

int rgb_match 
  ( long r, long g, long b, int maxval, 
    uint16_image_RGB_hist_vector cm, 
    int ncolors, 
    int mapmaxval
  );
  /* Search colormap for closest match in RGB metric. */
  
int yuv_match
  ( long r, long g, long b, int maxval, 
    uint16_image_RGB_hist_vector cm, 
    int ncolors, 
    int mapmaxval
  );
  /* Search colormap for closest match in YUV metric. */

#endif
