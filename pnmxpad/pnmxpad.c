#define PROG_NAME "pnmxpad"
#define PROG_DESC "add border to sides of a PBM/PGM/PPM file"
#define PROG_VERS "1.0"

/* Copyright © 2004 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2024-12-21 11:59:14 by stolfi */

#define PROG_HELP \
  PROG_NAME "  \\\n" \
  "    [ -white | -black | -color {COLOR} | -self ] \\\n" \
  "    [ -left {NUM} ] [ -right {NUM} ] [ -top {NUM} ] [ -bottom {NUM} ] \\\n" \
  "    [ {PNMFILE} ]\n" \
  "  where {COLOR} is\n" \
  "    " frgb_parse_color_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads a portable anymap as input, and" \
  " outputs a copy of the same, with extra rows and columns" \
  " around the border.\n" \
  "\n" \
  "  This program is similar to \"pnmpad\", with extra" \
  " features and different option syntax.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -white\n" \
  "  -black\n" \
  "  -color {COLOR}\n" \
  "    These options specify the value of the extra" \
  " pixels.  The {COLOR} argument is " frgb_parse_color_INFO "\n" \
  "\n" \
  "  -self\n" \
  "    This option specifies that each border pixel" \
  " should be equal to the nearest pixel of the" \
  " original image.\n" \
  "\n" \
  "  -left {NUM}\n" \
  "  -right {NUM}\n" \
  "  -top {NUM}\n" \
  "  -bottom {NUM}\n" \
  "    These options specify the width (in pixels)" \
  " of the extra border along each edge of the" \
  " image.  The default for each side is 0 (no extra border).\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmpad(1), ppmmake(1), pbmmake(1), pnmmargin(1), pnmpaste(1), pbm(5).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on jul/2004 by Jorge Stolfi, IC-UNICAMP.  \n" \
  "\n" \
  "  Based on {pnmpad.c} version 04/sep/1990 from" \
  " package {netpbm-1mar1994}, created by Jef Poskanzer (1989)" \
  " and Angus Duggan (1990).\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2004-07-31 Options \"-self\" and \"-color\" added by J. Stolfi, Unicamp.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2004 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <jsfile.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <argparser.h>

#define MAXV PNM_FILE_MAX_MAXVAL

/* COMMAND-LINE OPTIONS */

typedef struct rgb_pixel_t { uint16_t c[3]; } rgb_pixel_t; 

typedef struct options_t
  { char *fname; 
    bool_t self; 
    frgb_t bdcolor;
    int32_t left; 
    int32_t right;
    int32_t top;
    int32_t bot;
  } options_t;

/* PROTOTYPES */

int32_t main(int32_t argc, char **argv);
options_t *parse_options(int32_t argc, char **argv);
uint16_t *make_border_row(int32_t cols, int32_t chns, bool_t self, rgb_pixel_t bdcolor);
void fill_side_borders(uint16_t *smp, int32_t left, int32_t cols, int32_t right, int32_t chns, bool_t self, rgb_pixel_t bdcolor);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = parse_options(argc, argv);

    /* Open input file and get image/format paramters: */
    FILE *ifd = open_read(o->fname, FALSE);
    int32_t cols, rows, chns;
    uint16_t maxval;
    bool_t raw, bits;
    pnm_format_t format;
    pnm_read_header(ifd, &cols, &rows, &chns, &maxval, &raw, &bits, &format);
    
    if ((o->self) && ((cols == 0) || (rows == 0))) 
      { pnm_error("cannot self-pad an empty bitmap"); }
    
    /* Compute new dimensions: */
    int32_t newcols = cols + o->left + o->right;
    int32_t newrows = rows + o->top + o->bot;

    /* Convert border color from float to {0..maxval}: */
    rgb_pixel_t bdcolor;
    bool_t isMask = FALSE; /* Quantize as a smooth-valued image. */
    bdcolor.c[0] = pnm_quantize(o->bdcolor.c[0], maxval, isMask, PNM_NO_BADVAL);
    bdcolor.c[1] = pnm_quantize(o->bdcolor.c[1], maxval, isMask, PNM_NO_BADVAL);
    bdcolor.c[2] = pnm_quantize(o->bdcolor.c[2], maxval, isMask, PNM_NO_BADVAL);
    
    /* Allocate row buffer for new image: */
    uint16_t *imrow = uint16_image_alloc_pixel_row(newcols, chns);

    /* Allocate row buffer with given color: */
    uint16_t *bdrow = make_border_row(newcols, chns, o->self, bdcolor);

    pnm_write_header(stdout, newcols, newrows, maxval, format);

    for (uint32_t row = 0;  row < rows; row++) 
      { pnm_read_pixels(ifd, imrow + chns*o->left, cols, chns, maxval, raw, bits);
        fill_side_borders(imrow, o->left, cols, o->right, chns, o->self, bdcolor);

        /* Write top border before row 0: */
        if (row == 0)
          { for (uint32_t k = 0;  k < o->top; k++)
              { uint16_t *otrow = (o->self ? imrow : bdrow);
                pnm_write_pixels(stdout, otrow, newcols, chns, maxval, raw, bits);
              }
          }
          
        pnm_write_pixels(stdout, imrow, newcols, chns, maxval, raw, bits);

        /* Write bottom border after row {rows-1}: */
        if (row == rows-1)
          { for (uint32_t k = 0;  k < o->bot; k++)
              { uint16_t *otrow = (o->self ? imrow : bdrow);
                pnm_write_pixels(stdout, otrow, newcols, chns, maxval, raw, bits);
              }
          }
      }

    if (ifd != stdin) { fclose(ifd); }

    return 0;
  }

uint16_t *make_border_row(int32_t cols, int32_t chns, bool_t self, rgb_pixel_t bdcolor)
  { uint16_t *bdrow;
    if (self)
      { /* We do not need {bdrow}: */
        bdrow = NULL;
      }
    else
      { /* Make {bdrow} into a row of {o->bdcolor} pixels: */
        bdrow = uint16_image_alloc_pixel_row(cols, chns);
        int32_t col, c;
        uint16_t *p = bdrow;
        for (col = 0; col < cols; col++) 
          for (c = 0; c < chns; c++) 
            { (*p) = bdcolor.c[c]; p++; }
      }
    return bdrow;
  }        
 
void fill_side_borders(uint16_t *smp, int32_t left, int32_t cols, int32_t right, int32_t chns, bool_t self, rgb_pixel_t bdcolor)
  {
    uint16_t *p; /* Current pixel. */
    uint16_t *q; /* Source pixel. */
    /* Fill/replicate corner pixels into left/right borders: */
    q = smp + left*chns; /* Left edge. */
    int32_t col, c;
    for (col = 0, p = smp; col < left; col++) 
      for (c = 0; c < chns; c++) 
        { (*p) = (self ? (*(q+c)) : bdcolor.c[c]); p++; }
    q = smp + (left+cols-1)*chns; /* Right edge. */
    for (col = 0; col < right; col++) 
      for (c = 0; c < chns; c++) 
        { (*p) = (self ? (*(q+c)) : bdcolor.c[c]); p++; }
  }

options_t *parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    /* Parse keyword parameters: */
    if (argparser_keyword_present(pp, "-self"))
      { o->self = TRUE; o->bdcolor = frgb_Black; }
    else
      { if (argparser_keyword_present(pp, "-black")) 
          { o->self = FALSE; o->bdcolor = frgb_Black; }
        else if (argparser_keyword_present(pp, "-white")) 
          { o->self = FALSE; o->bdcolor = frgb_White; }
        else if (argparser_keyword_present(pp, "-color")) 
          { o->self = FALSE; o->bdcolor = frgb_parse_color(pp); }
        else
          { o->self = TRUE; o->bdcolor = frgb_Black; }
      }
      
    if (argparser_keyword_present(pp, "-left"))
      { o->left = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->left = 0; }
    
    if (argparser_keyword_present(pp, "-right"))
      { o->right = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->right = 0; }
    
    if (argparser_keyword_present(pp, "-bottom"))
      { o->bot = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->top = 0; }
    
    if (argparser_keyword_present(pp, "-top"))
      { o->top = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else
      { o->bot = 0; }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    if (argparser_next(pp) != NULL)
      { o->fname = argparser_get_next(pp); }
    else
      { o->fname = "-"; }

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
