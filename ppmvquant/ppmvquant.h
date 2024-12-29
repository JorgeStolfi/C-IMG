/* ppmvquant.h - options, types and prototypes for ppmvquant.c */
/* Last edited on 2024-12-26 19:51:13 by stolfi */

#ifndef ppmvquant_H
#define ppmvquant_H

#include <uint16_image_RGB_hist.h>
#include <uint16_image_RGB_table.h>
#include <uint16_color_tree.h>
#include <yuvhacks.h>
#include <uint16_image.h>
#include <uint16_image_match.h>

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

typedef struct box
  { int32_t ind;
    uint32_t colors;
    long sum;
  } box;
  
typedef box* box_vector;

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { bool_t floyd;
    char *map; 
    uint32_t maxColors; 
    uint16_image_match_proc_t *match;
    bool_t verbose;
    char *imgname;
  }options_t;

/* PROTOTYPES */

int32_t main(int32_t argc, char**argv);

options_t *parse_options(int32_t argc, char **argv);

uint32_t ppmvquant_choose_color
  ( int32_t r, int32_t g, int32_t b, 
    uint16_t maxval,
    uint16_image_RGB_hist_t *chv, 
    uint32_t nColors,
    uint16_t mapmaxval,
    uint16_image_match_proc_t *match,
    uint16_image_RGB_table_node_t **cht, 
    bool_t *addtohash
  );
  /* Finds best match to (r,g,b) in colormap chv, using/updating the 
    hash table if appropriate. */

void simple_quantize_uint16_image
  ( uint16_image_t *img, 
    uint16_image_RGB_hist_t *chv, 
    uint32_t newcolors,
    uint16_t mapmaxval,
    uint16_image_match_proc_t *match
  );

#endif
