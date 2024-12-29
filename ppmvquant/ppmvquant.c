#define PROG_NAME "ppmvquant"
#define PROG_DESC "quantize the colors in a PPM file with YUV or RGB metric"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-26 19:22:51 by stolfi */

/* Copyright © 1989, 1991 by Jef Poskanzer.
** See the copyright, authorship, and warranty notice at end of file.
*/

/* TO DO: !!! Update end fix -- syntax errors !!! */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [ -floyd {FLOYD} ] \\\n" \
  "  [ -yuv | -rgb ] \\\n" \
  "  { -maxColors {MAXCOLORS} | -map {MAPFILE} }\\\n" \
  "  [ {PPMFILE} ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program Reads a portable pixmap as input.  It" \
  " chooses up to {MAXCOLORS} colors to best represent the" \
  " image, maps the existing colors to the new ones, and" \
  " writes a portable pixmap as output.\n" \
  "\n" \
  "  Alternately, you can skip the color-choosing" \
  " step by specifying your own set of colors with" \
  " the \"-map\" argument.  The {MAPFILE} is just" \
  " a \".ppm\" file; the shape and size is" \
  " immaterial, all that matters is the DISTINCT" \
  " colors in it.   For instance, to quantize" \
  " down to the 8-color IBM TTL color set, you might use:\n" \
  "\n" \
  "      P3\n" \
  "      8 1\n" \
  "      255\n" \
  "        0   0   0\n" \
  "      255   0   0\n" \
  "        0 255   0\n" \
  "        0   0 255\n" \
  "      255 255   0\n" \
  "      255   0 255\n" \
  "        0 255 255\n" \
  "      255 255 255\n" \
  "\n" \
  "  Thus, if you want to quantize one pixmap to use the colors in" \
  " another one, just use the second one as the mapfile.   You don't" \
  " have to reduce it down to only one pixel of each color, just use it as is.\n" \
  "\n" \
  "  The quantization algorithm is Heckbert's /median cut/.  (See" \
  " Paul Heckbert (1982), \"Color Image Quantization for" \
  " Frame Buffer Display\"  SIGGRAPH 1982" \
  " Proceedings, page 297.) The nearest metch is" \
  " found using a visual metric (luminance first, then hue).\n" \
  "\n" \
  "  The \"-floyd\" option enables a Floyd-Steinberg error diffusion" \
  " step.  This algorithm gives vastly better results on images" \
  " where the unmodified quantization has banding or other" \
  " artifacts, especially when going to a small number of colors" \
  " such as the above IBM set.  However, it does take " \
  "substantially more CPU time, so the default is off.\n" \
  "\n" \
  "  This program is identical to Jef Poskanzer's" \
  " ppmquant(1), except for the option to use the" \
  " distance in YUV space when selecting the best match.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -floyd {FLOYD}\n" \
  "    This optional argument requests the Floyd-Steinberg" \
  " erro diffusion algorithm. The {FLOD} value may" \
  " be 'T', 'F', 0, or 1.  The default is \"-floyd F\".\n" \
  "\n" \
  "   -yuv\n" \
  "   -rgb\n" \
  "     These mutually exclusive optional arguments" \
  " specify the color space used to match colors.  The" \
  " default is \"-rgb\".\n" \
  "\n" \
  "   -maxColors {MAXCOLORS}\n" \
  "     Specifies the max number of distinct colors to" \
  " use in the output file.  It excludes the \"-map\" argument.\n" \
  "\n" \
  "   -map {MAPFILE}\n" \
  "     Specifies the name of a \".ppm\" file whose" \
  " distinct pixel colors are the color map to use" \
  " for the output file. It excludes the \"-maxColors\" argument.\n" \
  "  -blabber {AMOUNT}\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ppmquant(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in 1989-1991 by Jef Poskanzer\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  1996-05-18: YUV color metric added by J. Stolfi.\n" \
  "  1996-11-22: general rewrite to use JS libs by J. Stolfi.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 1989-1991 by Jef Poskanzer.\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_match.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>
#include <uint16_image_RGB_hist.h>
#include <uint16_image_RGB_table.h>
#include <uint16_image_RGB_medcut.h>
#include <uint16_image_quantize_floyd.h>
#include <argparser.h>
#include <affirm.h>

#include <ppmvquant.h>

uint16_image_RGB_hist_t choose_colormap 
  ( uint16_image_t *img,
    uint32_t maxColors,
    uint32_t *nColors_P
  );
  /*  Chooses a colormap with at most {maxColors} entries for {img}. 
    On exit, {*nColors_P} will be the number of colors actually used */
    
uint16_image_RGB_hist_t ppmvquant_collect_colors 
  ( uint16_image_t *map,
    uint32_t maxColors,
    uint32_t *nColors_P
  );    
  /* Collects all distinct colors in the image {map}. */
  
/* IMPLEMENTATIONS */

int main(int argc, char* argv[])
  { 
    options_t *o = parse_options(argc, argv);
    
    uint16_image_t *img = uint16_image_read_pnm_named(o->imgname, o->verbose);
    
    uint16_image_t *map;
    uint16_t mapmaxval;
    uint16_image_RGB_hist_t chv;
    uint32_t nColors; /* Actual number of colors in histogram {chv}. */
    if (o->map == NULL)
      { if (o->verbose) { fprintf(stderr, "using %d colors\n", o->maxColors); }
        chv = choose_colormap(img, o->maxColors, &nColors);
        mapmaxval = img->maxval;
      }
    else
      { if (o->verbose) { fprintf(stderr, "using given colormap\n"); }
        map = uint16_image_read_pnm_named(o->map, FALSE);
        mapmaxval = map->maxval;
        chv = ppmvquant_collect_colors(map, o->maxColors, &nColors);
      }
    
    set_yuv_matrix();
    if (o->floyd)
      { uint16_image_quantize_floyd(img, &chv, mapmaxval, o->match); }
    else
      { simple_quantize_uint16_image(img, &chv, mapmaxval, o->match); }

    uint16_image_write_pnm_named("-", img, 0, FALSE);
    return 0;
  }

uint16_image_RGB_hist_t choose_colormap (uint16_image_t *img, uint32_t maxColors, uint32_t *nColors_P)
  { demand(img->chns == 3, "image must have 3 channels");
    uint32_t nColors = maxColors;
    /* Attempt to make a histogram {chv0} of the colors, unclustered.
      If at first we don't succeed, lower maxval to increase color
      coherence and try again.  This will eventually terminate, with
      maxval at worst 15, since 32^3 is approximately MAXCOLORS. */
    uint16_image_RGB_hist_t chv0;
    while (TRUE)
      { chv0 = uint16_image_RGB_hist_build(img->smp, img->chns, img->cols, img->rows, &nColors);
        if (chv0.ne != 0) { break; }
        /* Reduce colors and retry: */
        uint16_t newmaxval = img->maxval / 2;
        demand(newmaxval > 0, "can't reduce {maxval} any further");
        int32_t row;
        uint32_t smp;
        uint16_t *pp;
        uint32_t spr = (uint32_t)(img->cols*img->chns);
        pnm_message( "too many colors, reducing maxval from %d to %d...", img->maxval, newmaxval);
        for (row = 0; row < img->rows; row++)
          for (smp = 0, pp = img->smp[row]; smp < spr; smp++)
            { double v = (((double)(*pp))*newmaxval)/img->maxval;
              uint32_t newsmp = (uint32_t)floor(v + 0.5);
              assert(newsmp <= newmaxval);
              (*pp) = (uint16_t)newsmp;
              pp++;
            }
        img->maxval = newmaxval;
      }
    (*nColors_P) = nColors;
    
    /* Now apply median-cut to histogram, making the new colormap. */
    uint16_image_RGB_hist_t chv1 = uint16_image_RGB_median_cut(&chv0, nColors, img->maxval, nColors_P);
    pnm_message("%d colors chosen out of %d", (*nColors_P), tcolors);

    free(chv0.e);
    return (chv1);
  }

uint16_image_RGB_hist_t ppmvquant_collect_colors (uint16_image_t *map, uint32_t *nColors_P)
  { if (map->cols == 0 || map->rows == 0) { pnm_error("null colormap??"); }
    uint16_image_RGB_hist_t chv = uint16_image_RGB_hist_build(map->smp, map->chns, map->cols, map->rows, nColors_P);
    if (chv.ne == 0) pnm_error("too many colors in colormap!");
    pnm_message("%d colors found in colormap", (*nColors_P));
    return(chv);
  }

void simple_quantize_uint16_image
  ( uint16_image_t *img, 
    uint16_image_RGB_hist_t *chv, 
    uint32_t maxColors,
    uint16_t mapmaxval,
    uint16_image_match_proc_t *match
  )
  { int32_t chns = img->chns; assert(chns == 3);
    int32_t cols = img->cols;
    int32_t rows = img->rows;
    bool_t addtohash = TRUE;
    uint16_image_RGB_table_node_t **cht = uint16_image_RGB_table_alloc();
    for (int32_t row = 0; row < rows; row++)
      { uint16_t *pp = img->smp[row];
        for (int32_t col = 0; col < cols; col++)
          { /* Choose replacement color* */
            int32_t sr = pp[0];
            int32_t sg = pp[1];
            int32_t sb = pp[2];
            ppm_pixel_t *qq = choose_color(
              sr, sg, sb, img->maxval,
              chv, maxColors, mapmaxval,
              match, cht, &addtohash
            );
            pp[0] = qq->c[0];
            pp[1] = qq->c[1];
            pp[2] = qq->c[2];
            pp += chns;
          }
      }
    img->maxval = mapmaxval;
  }


uint32_t ppmvquant_choose_color
  ( int32_t r, int32_t g, int32_t b,
    uint16_t maxval,
    uint16_image_color_hist_t *chv, 
    uint32_t *nColors,
    uint16_t mapmaxval,
    uint32_t maxColors,
    uint16_image_match_proc_t *match,
    uint16_color_table_node_t **cht
  )
  { uint32_t NH = (*nColors);
    /* Check whether the color is in the allowed range for search: */
    bool_t small = 
      (r >= 0) && (r <= uint16_image_MAX_MAXVAL) &&
      (g >= 0) && (g <= uint16_image_MAX_MAXVAL) &&
      (b >= 0) && (b <= uint16_image_MAX_MAXVAL);
    ppm_pixel_t p;
    int32_t ind = -1; /* Index into {chv} of ersatz color for {p}. */
    if (small)
      { /* Check whether we already looked up this color: */
        p = (ppm_pixel_t){{ (uint16_t)r, (uint16_t)g, (uint16_t)b }};;
        ind = uint16_color_table_lookup(cht, &p);
      }
    if (ind == -1)
      { /* Color is new, or is out of range due to error propagation. */
        /* Find the closest match: */
        ind = match(r, g, b, maxval, chv, NH, mapmaxval, cht);
        assert((ind >= 0) && (ind < chv->ne));
        if (small && (nColors < maxColors))
          if (uint16_color_table_add(cht, &p, ind) < 0)
            { pnm_message("no memory to grow hash table");
              (*addtohash) = 0;
            }
      }
    assert(ind >= 0);
    return(&(chv->e[ind].color));
  }


options_t *parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    if (argparser_keyword_present(pp, "-floyd"))
      { o->floyd = argparser_get_next_bool(pp); }
    else
      { o->floyd = TRUE; }
      
    if (argparser_keyword_present(pp, "-yuv"))
      { o->match = &uint16_image_match_yuv; }
    else if (argparser_keyword_present(pp, "-rgb"))
      { o->match = &uint16_image_match_rgb; }
    else
      { o->match = &uint16_image_match_rgb; }
      
    if (argparser_keyword_present(pp, "-maxColors"))
      { o->maxColors = (uint32_t)argparser_get_next_int(pp, 1, INT32_MAX);
        o->map = NULL;
      }
    else if (argparser_keyword_present(pp, "-map"))
      { o->map = argparser_get_next_non_keyword(pp);
        o->maxColors = 0; }
    else
      { argparser_error(pp, "either \"-maxColors\" or \"-map\" must be specified"); }
      
    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    if (argparser_next(pp) != NULL)
      { o->imgname = argparser_get_next(pp); }
    else
      { o->imgname = "-"; }

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
