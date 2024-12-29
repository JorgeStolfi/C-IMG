#define PROG_NAME "ppmxramp"
#define PROG_DESC "create a color ramp by interpolating three given pixels"
#define PROG_VERS "1.0"

/* Copyright © 2002 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2024-12-21 11:59:12 by stolfi */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -size {XSIZE} {YSIZE}  \\\n" \
  "    [ -maxval {MAXVAL} ] \\\n" \
  "    [ -channels {NC} ] \\\n" \
  "    {X[1]} {Y[1]} {COLOR[1]} \\\n" \
  "    {X[2]} {Y[2]} {COLOR[2]} \\\n" \
  "    {X[3]} {Y[3]} {COLOR[3]} \\\n" \
  "    > {OUTFILE} \\\n" \
  "  where each {COLOR[i]} is \\\n" \
  "    " frgb_parse_color_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program writes to standard output a PPM or PGM" \
  " image file containing a linear ramp (``gradient'') of" \
  " color values.  The ramp is such that the pixel in" \
  " column {X[i]} and row {Y[i]} has color {COLOR[i]}," \
  " for {i} in 1,2,3, and interpolates linearly between" \
  " those data points.\n" \
  "\n" \
  "  Pixel coordinates are measured from the upper left" \
  " corner of the image, with y increasing downwards.  The three" \
  " given colors apply to the center of the respective" \
  " pixels.  When generating grayscale output (with the" \
  " option \"-channels 1\"), the colors are coverted to" \
  " grayscale values by the formula {0.299*R + 0.587*G + 0.114*B}," \
  " as in {ppmtopgm}.\n" \
  "\n" \
  "  Each {COLOR[i]} is " frgb_parse_color_INFO "\n" \
  "\n" \
  "  The three {COLOR[i]} components are interpreted as RGB" \
  " coordinates in the [0 _ 1] range.\n" \
  "\n" \
  "  The output image file has linear encoding ({gamma == 1}).  Therefore," \
  " it will not produce a linear ramp of colors when displayed" \
  " on a typical monitor; and the hue may vary in an unexpected" \
  " manner over the image.  This ``feature'' can be corrected by piping the" \
  " output through {pnmgamma(1)}.\n" \
  "\n" \
  "  The ramp extends outside the triangle defined by the" \
  " three given pixels, and therefore the colors may saturate.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -size {XSIZE} {YSIZE}\n" \
  "    This mandatory argument defines the image dimensions, namely" \
  " {XSIZE} columns by {YSIZE} rows.\n" \
  "\n" \
  "  -maxval {MAXVAL}\n" \
  "    This optional argument specifies the maximum sample" \
  " value in the output image.  The default is the value of" \
  " {PNM_FILE_MAX_MAXVAL} in the interface {jspnm.h}.\n" \
  "\n" \
  "  -channels {NC}\n" \
  "    This optional argument specifies the number of channels" \
  " in the output image, which must be either 1 (grauscale) or" \
  " 3 (RGB color).  The default is 3.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmfield(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on aug/2002 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2008-07-07 by J. Stolfi, IC-UNICAMP.  Folded the manpage into the program.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2002 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <jspnm.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <uint16_image.h>
#include <argparser.h>

typedef struct rgb_pixel_t { uint16_t c[3]; } rgb_pixel_t; 

typedef struct options_t
  { uint16_t maxval;
    int32_t channels;
    double x[3];
    double y[3];
    frgb_t v[3];
    int32_t rows;
    int32_t cols;
  } options_t;

options_t *get_options(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    int32_t chns = o->channels;
    bool_t raw, bits;
    pnm_format_t format;
    pnm_choose_output_format(o->maxval, chns, FALSE, &format, &raw, &bits);
    pnm_write_header(stdout, o->cols, o->rows, o->maxval, format);
    uint16_t *smprow = uint16_image_alloc_pixel_row(o->cols, chns);

    /* Compute the Cartesian-to-barycentric coordinate transform matrix {M}: */
    double M11, M12, M21, M22;
    { double a11 = o->x[1] - o->x[0];
      double a12 = o->x[2] - o->x[0];
      double a21 = o->y[1] - o->y[0];
      double a22 = o->y[2] - o->y[0];
      double deta = a11*a22-a21*a12;
      if (deta == 0)
        { fprintf(stderr, "%s: ** reference points are collinear\n", PROG_NAME);
          exit(1);
        }
      M11 = +a22/deta;
      M12 = -a12/deta;
      M21 = -a21/deta;
      M22 = +a11/deta;
    }

    auto void cart_to_bary(int32_t col, int32_t row, double alpha[]);
      /* Computes the barycentric coordinates {alpha[0..2]} of pixel {(col,row)}: */

    void cart_to_bary(int32_t col, int32_t row, double alpha[])
      { double bb1  = (double)col - o->x[0];
        double bb2  = (double)row - o->y[0];
        alpha[1] = bb1*M11+bb2*M12;
        alpha[2] = bb1*M21+bb2*M22;
        alpha[0] = 1.0 - alpha[1] - alpha[2];
      }

    double alpha[chns]; /* Barycentric coordinates of pixel rel. to ref pixels. */

    if (chns == 1)
      { double y0 = frgb_get_Y_pbm(&(o->v[0]));
        double y1 = frgb_get_Y_pbm(&(o->v[1]));
        double y2 = frgb_get_Y_pbm(&(o->v[2]));

        for (int32_t row = 0; row < o->rows; ++row)
          { uint16_t *sp = smprow;
            for (int32_t col = 0; col < o->cols; ++col)
              { cart_to_bary(col, row, alpha);
                double Y = alpha[0]*y0 + alpha[1]*y1 + alpha[2]*y2;
                if (Y < 0) { Y = 0; }
                if (Y > 1) { Y = 1; }
                int32_t q = (int32_t)floor(Y*o->maxval + 0.5);
                (*sp) = (uint16_t)q; sp++;
              }
            pnm_write_pixels(stdout, smprow, o->cols, chns, o->maxval, raw, bits);
          }
      }
    else if (chns == 3)
      { frgb_t *v0 = &(o->v[0]);
        frgb_t *v1 = &(o->v[1]);
        frgb_t *v2 = &(o->v[2]);

        /* Compute the mean {B} of the corner colors: */
        frgb_t B;
        for (uint32_t chn = 0;  chn < chns; chn++)
          { B.c[chn] = (float)(((double)v0->c[chn] + (double)v1->c[chn] + (double)v2->c[chn])/3.0); }

        for (int32_t row = 0; row < o->rows; ++row)
          { uint16_t *sp = smprow;
            for (int32_t col = 0; col < o->cols; ++col)
              { cart_to_bary(col, row, alpha);
                frgb_t S;
                for (uint32_t chn = 0;  chn < chns; chn++)
                  { S.c[chn] = (float)(alpha[0]*v0->c[chn] + alpha[1]*v1->c[chn] + alpha[2]*v2->c[chn]); }
                frgb_clip_rgb_towards(&S, &B);
                for (uint32_t chn = 0;  chn < chns; chn++)
                  { int32_t q = (int32_t)floor(S.c[chn]*o->maxval + 0.5);
                    (*sp) = (uint16_t)q; sp++;
                  }
              }
            pnm_write_pixels(stdout, smprow, o->cols, chns, o->maxval, raw, bits);
          }
      }

    fflush(stdout);

    return 0;
  }

#define MAX_IMAGE_SIZE (16*1024)
  /* To avoid humongous allocs. */

options_t *get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    
    /* Parse the keyword arguments: */
    
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (uint16_t)argparser_get_next_int(pp, 1, PNM_FILE_MAX_MAXVAL); }
    else
      { o->maxval = PNM_FILE_MAX_MAXVAL; }
      
    if (argparser_keyword_present(pp, "-channels"))
      { o->channels = (int32_t)argparser_get_next_int(pp, 1, 3);
        if ((o->channels != 1) && (o->channels != 3))
          { argparser_error(pp, "the number of channels must be 1 or 3"); }
      }
    else
      { o->channels = 3; }
   
    argparser_get_keyword(pp, "-size");
    o->cols = (int32_t)argparser_get_next_int(pp, 0, MAX_IMAGE_SIZE);
    o->rows = (int32_t)argparser_get_next_int(pp, 0, MAX_IMAGE_SIZE);

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
 
    for (uint32_t i = 0;  i < 3; i++)
      { o->x[i] = argparser_get_next_double(pp, -1.0e+100, +1.0e+100);
        o->y[i] = argparser_get_next_double(pp, -1.0e+100, +1.0e+100);
        o->v[i] = frgb_parse_color(pp);
      }

    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }
