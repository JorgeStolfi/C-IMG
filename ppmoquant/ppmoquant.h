#ifndef ppmoquant_H
#define ppmoquant_H

/* ppmoquant.h - options, types and prototypes for ppmoquant.c */
/* Last edited on 2017-06-20 20:43:05 by stolfilocal */

#include <jspnm.h>
#include <uint16_image.h>
/* #include <jsppmcmap.h> */
#include <uint16_image_RGB_table.h>
#include <argparser.h>

#define MAXCOLORS 32767

/* COMMAND-LINE OPTIONS */

typedef int matchfn_t
  ( long br, long bg, long bb, 
    long wr, long wg, long wb, 
    int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval
  );
  /* Type of a function that finds the best (opaque) match to the color
    interval [r1_r2]x[g1_g2]x[b1_b2] in the colormap cm. Note that 
    the interval may be somewhat outside the color cube. */

typedef struct options_t
  { bool_t floyd; 
    matchfn_t *match;
    ppm_pixel_t transpcolor;
    int transpmaxval;
    char *cmapname;
    char *bimgname; 
    char *wimgname;
    bool_t verbose;
  } options_t;

/* Prototypes */

int main(int argc, char* argv[]);

options_t *parse_options(int argc, char **argv);

ppm_pixel_t my_parsecolor(argparser_t *pp);
  /* Parses three integer arguments from the commad line,
    returns them as the three samples of a {ppm_pixel_t}. */

ppm_pixel_t* find_color
  ( ppm_pixel_t p,
    uint16_image_RGB_hist_vector cm, 
    int ncolors
  );
  /* Looks for an exact match to p in cm[0..ncolors-1].
    Returns pointer to pixel if found, NULL if not found. */

uint16_image_RGB_hist_vector collect_colors(uint16_image_t *map, int *ncolorsp);
  /* Builds from {map->samples} a color histogram {uint16_image_RGB_hist_vector}. */

void uint16_image_quantize_plains
  ( uint16_image_t *bimg, 
    uint16_image_t *wimg, 
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    ppm_pixel_t transpcolor,
    int mapmaxval,
    matchfn_t *match
  );
  /* Replaces each pixel in {bimg} by a pixel from the given colormap,
    or the "transpcolor" (transparent) pixel, by comparing the two colors. */

void uint16_image_quantize_floyds
  ( uint16_image_t *bimg, 
    uint16_image_t *wimg, 
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    ppm_pixel_t transpcolor,
    int mapmaxval,
    matchfn_t *match
  );
  /* Replaces each pixel in {bimg} by a pixel from the given colormap,
    or the {transpcolor} (transparent) pixel, by comparing the two colors.
    Uses Floyd-Steinberg error propagation. */

typedef struct
  { long* thisr;
    long* nextr;
    long* thisg;
    long* nextg;
    long* thisb;
    long* nextb;
  } fs_errors;

void floyd_alloc_error_vectors(fs_errors *ep, int cols);

void floyd_propagate_error_opaque
  ( long br, long bg, long bb,
    long wr, long wg, long wb,
    ppm_pixel_t *newp, double coef,
    fs_errors *bep, fs_errors *wep,
    int col,
    int fs_direction
  );
  /* Propagates Floyd-Steinberg error terms for opaque pixel choice.
    The old pixels are br, bg, bb, wr, wg, wb.
    The new pixel is (*newp).  
    Assumes coef = FS_SCALE*(oldmaxval/newmaxval). */
  
void floyd_propagate_error_transp
  ( long br, long bg, long bb,
    long wr, long wg, long wb,
    int maxval,
    fs_errors *bep, fs_errors *wep,
    int col,
    int fs_direction
  );
  /* Propagates Floyd-Steinberg error terms for transparent pixel choice.
    The old pixels are br, bg, bb, wr, wg, wb. */
  
void floyd_swap_error_vectors(fs_errors *ep);

ppm_pixel_t* choose_color
  ( long r0, long g0, long b0,
    long r1, long g1, long b1,
    int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval,
    matchfn_t match
  );
  /* Finds best match to the color interval (r0,g0,b0)--(r1,g1,b1) 
    in colormap cm plus the "transparent" color.  Returns NULL
    if the best match is "transparent". */
  
int rgb_match_range
  ( long br, long bg, long bb,
    long wr, long wg, long wb,
    int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval
  );
  /* Search colormap for closest match to given color range. */

#endif
