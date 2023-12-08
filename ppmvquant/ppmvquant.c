#define PROG_NAME "ppmvquant"
#define PROG_DESC "quantize the colors in a PPM file with YUV or RGB metric"
#define PROG_VERS "1.0"

/* Copyright © 1989, 1991 by Jef Poskanzer.
** See the copyright, authorship, and warranty notice at end of file.
** Last edited on 2017-06-20 20:43:35 by stolfilocal
*/

/* TO DO: !!! Update end fix -- syntax errors !!! */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [ -fs | -nofs ] [ -yuv | -rgb ] \\\n" \
  "  { NCOLORS | -map MAPFILE } \\\n" \
  "  [ PPMFILE ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program bla bla bla" \
  " bla bla bla bla {X+Y} bla bla" \
  " bla {INFILE} bla \"foobar.ppm\" bla bla bla\n" \
  "\n" \
  "  Beware that bla bla bla BLEBBLE BLOB bla" \
  " bla bla bla bla.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -blabber {AMOUNT}\n" \
  "    Blabbers for that {AMOUNT}. May also bla" \
  " bla bla bla bla bla bla bla bla bla bla bla bla" \
  " bla bla bla bla bla bla bla.\n" \
  "\n" \
  "  This program is similar to \"ppmquant\" but uses YUV color metric as an option.\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ppmquant(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 1989-1991 by Jef Poskanzer.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  22/nov/1996: general rewrite to use JS libs by J. Stolfi.\n" \
  "\n" \
  "  18/may/1996: YUV color metric added by J. Stolfi.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 1989-1991 by Jef Poskanzer.\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS


#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_RGB_hist.h>
#include <uint16_image_RGB_table.h>
#include <uint16_image_RGB_medcut.h>
#include <argparser.h>
#include <affirm.h>

#include <ppmvquant.h>

/* IMPLEMENTATIONS */

int main(int argc, char* argv[])
  { 
    options_t *o = parse_options(argc, argv);
    
    uint16_image_RGB_hist_vector cm;

    uint16_image_t *img = uint16_image_read_pnm_named(o->imgname, o->verbose);
    
    uint16_image_t *map;
    int mapmaxval;
    if (o->mapname == NULL)
      { if (o->verbose) { fprintf(stderr, "using %d colors\n", o->ncolors); }
        cm = choose_colormap(img, &(o->ncolors));
        mapmaxval = img->maxval; }
    else
      { if (o->verbose) { fprintf(stderr, "using colormap\n"); }
        map = uint16_image_read_pnm_named(o->mapname, FALSE);
        mapmaxval = map->maxval;
        cm = collect_colors(map, &(o->ncolors));
      }
    
    set_yuv_matrix();
    if (o->floyd)
      { uint16_image_quantize_floyd(img, cm, o->ncolors, mapmaxval, o->match); }
    else
      { simple_quantize_uint16_image(img, cm, o->ncolors, mapmaxval, o->match); }

    uint16_image_write_pnm_named("-", img, 0, FALSE);
    return 0;
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
    
    if (argparser_keyword_present(pp, "-fs") || argparser_keyword_present(pp, "-floyd"))
      { o->floyd = TRUE; }
    else if (argparser_keyword_present(pp, "-nofs") || argparser_keyword_present(pp, "-nofloyd"))
      { o->floyd = FALSE; }
    else
      { o->floyd = TRUE; }
      
    if (argparser_keyword_present(pp, "-yuv"))
      { o->match = yuv_match; }
    else if (argparser_keyword_present(pp, "-rgb"))
      { o->match = rgb_match; }
    else
      { o->match = rgb_match; }
      
    if (argparser_keyword_present(pp, "-map"))
      { o->mapname = argparser_get_next(pp); }
    else
      { o->mapname = NULL; }

    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    if (o->mapname == NULL)
      { if (argparser_next(pp) == NULL) 
          { argparser_error(pp, "missing number of colors"); }
        o->ncolors = argparser_get_next_int(pp, 1, INT_MAX);
      }
    else
      { o->ncolors = 0; }

    if (argparser_next(pp) != NULL)
      { o->imgname = argparser_get_next(pp); }
    else
      { o->imgname = "-"; }

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

uint16_image_RGB_hist_vector choose_colormap (uint16_image_t *img, int *ncolorsp)
  /*
    Chooses a colormap with at most *ncolors entries for *img.
    On exit *ncolors may have been reduced.
  */
  { demand(img->chns == 3, "image must have 3 channels");
    int tcolors;
    /* Attempt to make a histogram {ch} of the colors, unclustered.
      If at first we don't succeed, lower maxval to increase color
      coherence and try again.  This will eventually terminate, with
      maxval at worst 15, since 32^3 is approximately MAXCOLORS. */
    uint16_image_RGB_hist_vector ch;
    while (TRUE)
      { ch = uint16_image_RGB_hist_build(img->smp, img->chns, img->cols, img->rows, MAXCOLORS, &tcolors);
        if (ch != NULL) { break; }
        /* Reduce colors and retry: */
        int newmaxval = img->maxval / 2;
        int row, smp;
        uint16_t *pp;
        int spr = img->cols * img->chns;
        pnm_message
          ( "too many colors, reducing maxval from %d to %d...", img->maxval, newmaxval);
        for (row = 0; row < img->rows; row++)
          for (smp = 0, pp = img->smp[row]; smp < spr; smp++)
            { double v = (((double)(*pp))*newmaxval)/img->maxval;
              (*pp) = (int)(v + 0.5);
              pp++;
            }
        img->maxval = newmaxval;
      }
    
    /* Now apply median-cut to histogram, making the new colormap. */
    uint16_image_RGB_hist_vector cm = uint16_image_RGB_median_cut(ch, tcolors, img->maxval, ncolorsp);
    pnm_message("%d colors chosen out of %d", (*ncolorsp), tcolors);

    uint16_image_RGB_hist_free(ch);
    return (cm);
  }

uint16_image_RGB_hist_vector collect_colors (uint16_image_t *map, int *ncolorsp)
  { if (map->cols == 0 || map->rows == 0) { pnm_error("null colormap??"); }
    uint16_image_RGB_hist_vector cm = uint16_image_RGB_hist_build(map->smp, map->chns, map->cols, map->rows, MAXCOLORS, ncolorsp);
    if (cm == NULL) pnm_error("too many colors in colormap!");
    pnm_message("%d colors found in colormap", (*ncolorsp));
    return(cm);
  }

void uint16_image_quantize_floyd
  ( uint16_image_t *img, 
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval,
    matchfn_t *match
  )
  { int chns = img->chns; assert(chns == 3);
    int cols = img->cols;
    int rows = img->rows;
    
#define FS_SCALE 1024
#define FS_WTA 7
#define FS_WTB 3
#define FS_WTC 5
#define FS_WTD 1
    long* thisrerr = (long*)notnull(malloc((img->cols + 2)*sizeof(long)), "no mem");
    long* nextrerr = (long*)notnull(malloc((img->cols + 2)*sizeof(long)), "no mem");
    long* thisgerr = (long*)notnull(malloc((img->cols + 2)*sizeof(long)), "no mem");
    long* nextgerr = (long*)notnull(malloc((img->cols + 2)*sizeof(long)), "no mem");
    long* thisberr = (long*)notnull(malloc((img->cols + 2)*sizeof(long)), "no mem");
    long* nextberr = (long*)notnull(malloc((img->cols + 2)*sizeof(long)), "no mem");
    long* temperr;
    int row, col;
    int maxval = img->maxval;
    long maxcor = 2*maxval;
    long mincor = -1*maxval;
    int fs_direction = 1;
    int addtohash = 1;
    uint16_image_RGB_table cht = uint16_image_RGB_table_alloc();
    long sr, sg, sb, new, err;
    float coef = FS_SCALE*((float)maxval)/((float)mapmaxval);

    /* Initialize Floyd-Steinberg error vectors with randoms in [-1..+1]. */
    srand((int) 46157);
    for (col = 0; col < cols + 2; col++)
      { thisrerr[col] = (rand() % (FS_SCALE * 2)) - FS_SCALE;
        thisgerr[col] = (rand() % (FS_SCALE * 2)) - FS_SCALE;
        thisberr[col] = (rand() % (FS_SCALE * 2)) - FS_SCALE;
      }

    /* Scan the image rows: */
    for (row = 0; row < rows; row++)
      { /* Clear the {nextXerr} arrays: */
        for (col = 0; col < cols + 2; col++)
          { nextrerr[col] = nextgerr[col] = nextberr[col] = 0; }
        /* Define the initial and final column */
        int limitcol;
        if (fs_direction)
          { col = 0; limitcol = cols; }
        else
          { col = cols - 1; limitcol = -1; }
        /* Scan columns and propagate errors */
        uint16_t *pp = &(img->smp[row][col*chns]);
        do
          { /* Use Floyd-Steinberg errors to adjust actual color. */
            sr = pp[0] + thisrerr[col + 1] / FS_SCALE;
            sg = pp[1] + thisgerr[col + 1] / FS_SCALE;
            sb = pp[2] + thisberr[col + 1] / FS_SCALE;
            if (sr < mincor) { sr = mincor; } if (sr > maxcor) { sr = maxcor; }
            if (sg < mincor) { sg = mincor; } if (sg > maxcor) { sg = maxcor; }
            if (sb < mincor) { sb = mincor; } if (sb > maxcor) { sb = maxcor; }
            /* Choose replacement color* */
            ppm_pixel_t *qq = choose_color(
              sr, sg, sb, maxval, 
              cm, ncolors, mapmaxval,
              match, cht, &addtohash
            );
            /* Replace pixel in image: */
            pp[0] = qq->c[0];
            pp[1] = qq->c[1];
            pp[2] = qq->c[2];
            /* Propagate Floyd-Steinberg error terms. */
            if (fs_direction)
              { new = (long)(coef*qq->c[0]+0.5);
                err = (sr * FS_SCALE - new);
                thisrerr[col + 2] += (err * FS_WTA) / 16;
                nextrerr[col    ] += (err * FS_WTB) / 16;
                nextrerr[col + 1] += (err * FS_WTC) / 16;
                nextrerr[col + 2] += (err * FS_WTD) / 16;
                new = (long)(coef*qq->c[1]+0.5);
                err = (sg * FS_SCALE - new);
                thisgerr[col + 2] += (err * FS_WTA) / 16;
                nextgerr[col    ] += (err * FS_WTB) / 16;
                nextgerr[col + 1] += (err * FS_WTC) / 16;
                nextgerr[col + 2] += (err * FS_WTD) / 16;
                new = (long)(coef*qq->c[2]+0.5);
                err = (sb * FS_SCALE - new);
                thisberr[col + 2] += (err * FS_WTA) / 16;
                nextberr[col    ] += (err * FS_WTB) / 16;
                nextberr[col + 1] += (err * FS_WTC) / 16;
                nextberr[col + 2] += (err * FS_WTD) / 16;
                col++; pp += chns;
              }
            else
              { new = (long)(coef*qq->c[0]+0.5);
                err = (sr * FS_SCALE - new);
                thisrerr[col    ] += (err * FS_WTA) / 16;
                nextrerr[col + 2] += (err * FS_WTB) / 16;
                nextrerr[col + 1] += (err * FS_WTC) / 16;
                nextrerr[col    ] += (err * FS_WTD) / 16;
                new = (long)(coef*qq->c[1]+0.5);
                err = (sg * FS_SCALE - new);
                thisgerr[col    ] += (err * FS_WTA) / 16;
                nextgerr[col + 2] += (err * FS_WTB) / 16;
                nextgerr[col + 1] += (err * FS_WTC) / 16;
                nextgerr[col    ] += (err * FS_WTD) / 16;
                new = (long)(coef*qq->c[2]+0.5);
                err = (sb * FS_SCALE - new);
                thisberr[col    ] += (err * FS_WTA) / 16;
                nextberr[col + 2] += (err * FS_WTB) / 16;
                nextberr[col + 1] += (err * FS_WTC) / 16;
                nextberr[col    ] += (err * FS_WTD) / 16;
                col--; pp -= chns;                
              }
          }
        while (col != limitcol);

        temperr = thisrerr;
        thisrerr = nextrerr;
        nextrerr = temperr;
        temperr = thisgerr;
        thisgerr = nextgerr;
        nextgerr = temperr;
        temperr = thisberr;
        thisberr = nextberr;
        nextberr = temperr;
        fs_direction = ! fs_direction;
      }
    img->maxval = mapmaxval;
  }

void simple_quantize_uint16_image
  ( uint16_image_t *img, 
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval,
    matchfn_t *match
  )
  { int chns = img->chns; assert(chns == 3);
    int cols = img->cols;
    int rows = img->rows;
    int row;
    long sr, sg, sb;
    int addtohash = 1;
    uint16_image_RGB_table cht = uint16_image_RGB_table_alloc();
    for (row = 0; row < rows; row++)
      { int col = 0;
        int limitcol = cols;
        uint16_t *pp = img->smp[row];
        do
          { /* Choose replacement color* */
            sr = pp[0];
            sg = pp[1];
            sb = pp[2];
            ppm_pixel_t *qq = choose_color(
              sr, sg, sb, img->maxval,
              cm, ncolors, mapmaxval,
              match, cht, &addtohash
            );
            pp[0] = qq->c[0];
            pp[1] = qq->c[1];
            pp[2] = qq->c[2];
            col++; pp += chns;
          }
        while (col != limitcol);
      }
    img->maxval = mapmaxval;
  }

ppm_pixel_t *choose_color
  ( long r, long g, long b, int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval,
    matchfn_t *match,
    uint16_image_RGB_table cht, 
    int *addtohash
  )
  { int ind = -1;
    int small = 0;
    ppm_pixel_t p;
    if ((r >= 0) && (r <= PNM_FILE_MAX_MAXVAL) &&
        (g >= 0) && (g <= PNM_FILE_MAX_MAXVAL) &&
        (b >= 0) && (b <= PNM_FILE_MAX_MAXVAL))
      { /* Check hash table to see if we have already matched this color. */
        small = 1;
        p = (ppm_pixel_t){{ r, g, b }};;
        ind = uint16_image_RGB_table_lookup(cht, &p);
      }
    if (ind == -1)
      { ind = match(r, g, b, maxval, cm, ncolors, mapmaxval);
        if (small && (*addtohash))
          if (uint16_image_RGB_table_add(cht, &p, ind) < 0)
            { pnm_message("no memory to grow hash table");
              (*addtohash) = 0;
            }
      }
    assert(ind >= 0);
    return(&(cm[ind].color));
  }

int rgb_match
  ( long r, long g, long b, int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval
  )
  { register int i;
    register int ind = -1;
    register float y1, u1, v1, s, y2, u2, v2;
    register float dist, newdist;
    s = 1.0/(float)maxval;
    y1 = s*r;
    u1 = s*g;
    v1 = s*b;
    dist = 1.0e20;
    s = 1.0/(float)mapmaxval;
    for (i = 0; ((dist > 0) && (i < ncolors)); i++)
      { ppm_pixel_t *qq = &(cm[i].color);
        y2 = s * qq->c[0];
        u2 = s * qq->c[1];
        v2 = s * qq->c[2];
        newdist = (y1-y2)*(y1-y2) + (u1-u2)*(u1-u2) + (v1-v2)*(v1-v2);
        if (newdist < dist) { ind = i; dist = newdist; }
      }
    return (ind);
  }

int yuv_match 
  ( long r, long g, long b, int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval
  )
  { register int i;
    register int ind = -1;
    int ri, gi, bi;
    float y1, u1, v1, y2, u2, v2;
    register float dist, newdist;
    rgb_to_yuv(r, g, b, maxval, &y1, &u1, &v1);
    dist = 1.0e20;
    for (i = 0; ((dist > 0) && (i < ncolors)); i++)
      { ppm_pixel_t *qq = &(cm[i].color);
        ri = qq->c[0];
        gi = qq->c[1];
        bi = qq->c[2];
        rgb_to_yuv(ri, gi, bi, mapmaxval, &y2, &u2, &v2);
        newdist = (y1-y2)*(y1-y2) + (u1-u2)*(u1-u2) + (v1-v2)*(v1-v2);
        if (newdist < dist) { ind = i; dist = newdist; }
      }
    return (ind);
  }
