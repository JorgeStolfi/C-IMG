#define PROG_NAME "make_bullseye"
#define PROG_DESC "creates a grayscale image with a circular wave pattern"
#define PROG_VERS "1.0"

/* Last edited on 2023-03-03 15:39:26 by stolfi */

/* Copyright © 2007 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -size {NX} {NY} \\\n" \
  "    > OUTFILE.pgm"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Writes to stdout a portable graymap (PGM) image file with a" \
  " circular wave pattern.  The wavelength is inversely proportional" \
  " to the distance from the center.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -size {NX} {NY}\n" \
  "    This mandatory argument specifies number of columns and number of rows of the image.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgm(5), and anything else you should see.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created around  Aug 16  2007-08-16 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-03-03 J.Stolfi: added documentation, option parsing, {stdint}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2007 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <argparser.h>
#include <jspnm.h>

typedef struct options_t 
  { int32_t size_NX;
    int32_t size_NY;
  } options_t;

int32_t main(int32_t argc, char **argv);
options_t *get_options(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);

    int32_t NX = o->size_NX;
    int32_t NY = o->size_NY;
    
    uint16_t maxval = PNM_FILE_MAX_MAXVAL;

    /* Write PGM (ascii) header: */
    fprintf(stdout,"P2\n");
    fprintf(stdout,"%d %d\n", NX, NY);
    fprintf(stdout,"%u\n", maxval);

    double S = hypot(NX, NY)/2;
    int32_t NH = 5; /* Pixel subsampling order: */
    for(int32_t i = 0; i < NY; i++)
      for(int32_t j = 0; j < NX; j++)
        { /* Compute pixel in row {i} and column {j}, antialiased: */ 
          double sum_fw = 0;
          double sum_w = 0;
          double cx = (j + 0.5 - 0.5*NX);
          double cy = (i + 0.5 - 0.5*NY);

          for (int32_t di = -NH+1; di <= +NH-1; di++) 
            for (int32_t dj = -NH+1; dj <= +NH-1; dj++) 
              { 
                double dx = ((double)dj)/NH;
                double wx = 0.5*(1 + cos(M_PI*dx));
                double zx = cx + dx;

                double dy = ((double)di)/NH;
                double wy = 0.5*(1 + cos(M_PI*dy));
                double zy = cy + dy;

                double r2 = zx*zx + zy*zy;
                double w = wx*wy;

                double f = sin(r2/S);
                sum_fw += f*w;
                sum_w += w;
              }
          double avg_f = sum_fw/sum_w;
          assert(fabs(avg_f) <= 1.0);
          double pval = 0.5*(1 + avg_f); /* Pixel value in {[0_1]}. */
          double pscale = 0.4999999*(maxval - 1);
          int32_t dval = (1 + (int32_t)floor(pscale*pval + 0.5));
          assert((dval >= 1) && (dval <= maxval));
          fprintf(stdout," %u\n", (uint16_t)dval);
        }
    fflush(stdout);
    return 0;
  }

options_t *get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 

    argparser_get_keyword(pp, "-size");
    o->size_NX = (int32_t)argparser_get_next_int(pp, 1, PNM_FILE_MAX_MAXVAL);
    o->size_NY = (int32_t)argparser_get_next_int(pp, 1, PNM_FILE_MAX_MAXVAL);

    argparser_skip_parsed(pp);
    
    argparser_finish(pp);

    return o;
  }

