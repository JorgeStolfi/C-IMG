#define PROG_NAME "optdither"
#define PROG_DESC "generates ordered dither matrices by spectrum optimization"
#define PROG_VERS "2.0"

/* Copyright © 2007 by the State University of Campinas (UNICAMP). */
/* See the authorship, rights and warranty notices in the PROG_INFO below. */
/* Last edited on 2023-03-18 09:13:51 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -limit {NITER} \\\n" \
  "    -blabber {AMOUNT} \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {OUTFILE}"

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
  "  -niter {NITER}\n" \
  "    This mandatory argument specifies the number of iterations.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ppmquant(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created by Konstantin Glasov in 1993 at IC-UNICAMP.\n" \
  "  Rewritten 2007-08-15 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  1993-11-30 Program {iterations.c} handed in by K.Glasov at the end of his stay.\n" \
  "  2007-08-15 Rewritten by J.Stolfi to use {argparser.h} and standard help/info.\n" \
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
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <gauss_table.h>
#include <uint16_image.h>
#include <uint16_image_read_gen.h>
#include <float_image_write_pnm.h>
#include <float_image.h>
#include <image_file_format.h>

#include <odt_tools.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { int32_t niter;
    bool_t op1;     /* A boolean argument. */
    char* op2;      /* A string-valued argument. */
  } options_t;

/* INTERNAL PROTOTYPES */

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int32_t main(int32_t argc,char** argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);
 
    /* Read the initial dither matrix {xo}: */
    char *fname = o->initDither;
    uint32_t maxval;
    float_image_t *xo = odt_read_dither_matrix(fname, &maxval);
   
    int32_t NC, NX, NY;
    float_image_get_size(xo, &NC, &NX, &NY);
    assert(NC == 1);
    
    /* Create the filter tables: */
    double fmid = o->filterFreq;
    bool_t norm = FALSE;  /* Don't normalize to unit sum. */
    bool_t folded = TRUE; /* Fold the gaussian  modulo {n}. */
    double *wtx = gauss_table_make(NX, 0.0, fmid, norm, folded);
    double *wty = gauss_table_make(NY, 0.0, fmid, norm, folded);

    /* Create the auxiliary images: */
    float_image_t *yo = float_image_new(NC, NX, NY);
    float_image_t *wo = float_image_new(NC, NX, NY);
    
    /* Optimization loop: */
    int32_t max_iter = o->maxIter;
    int32_t iter = 0;
    while (iter < max_iter)
      { iter++;

        odt_one_step(xo, xf, yf, yo, wtx, wty);
        if ((iter <= 10) || ((iter % 10) == 1))
          { odt_show_iteration(o->outPrefix, iter, xo, xf, yf, yo); }

        odt_compare_matrices(xo, yo);

        /* Swap the arrays: */
        { float_image_t *ro = xo; xo = yo; yo = ro; }
      }
    
    odt_show_matrix(o->outPrefix, -1, "_final_xo", xo);
   
    return 0;
  }
  
void odt_one_step
  ( float_image_t *xo, 
    float_image_t *xf, 
    float_image_t *yf, 
    float_image_t *yo,
    double wtx[],
    double wty[]
  )
  {
    /* Compute the Fourier transform of {x}: */
    float_image_hartley_transform(xo, xf);
    
    /* Reduce the low-freq components: */
    odt_fft_filter(xf, yf, wtx, wty);
    
    /* Invert the Fourier transform: */
    float_image_hartley_transform(yf, yo);
    
    /* Round off to a permutation: */
    odt_fft_permutize(yo);
    
  }

void odt_show_iteration
  ( char *prefix, 
    int32_t iter,
    float_image_t *xo, 
    float_image_t *xf, 
    float_image_t *yf, 
    float_image_t *yo
  )
  {
    odt_show_matrix(prefix, iter, "xo", xo);
    odt_show_ft(prefix, iter, "xf", xf);
    odt_show_ft(prefix, iter, "yf", yf);
    odt_show_matrix(prefix, iter, "yo", yo);
  }
  
void odt_show_matrix(char *prefix, iny32_t iter, char *tag, float_image_t *xo)
  {
    /* Convert to integer image: */
    assert(NC == 1);
    bool_t isMask = FALSE,
    uint32_t maxval = NX*NY-1;
    assert(maxval <= uint16_image_MAX_SAMPLE);
    bool_t yup = TRUE;
    bool_t verbose_conv = FALSE;
    uint16_image_t *xp = float_image_to_uint16_image
      ( xo, isMask, 1, NULL, NULL, NULL, maxval, yup, verbose_conv );

    /* Write to disk: */
    char *fname = NULL;
    if (iter >= 0)
      { asprintf(&fname, "%s_%04d_%s.pgm", prefix, iter, tag); }
    else
      { asprintf(&fname, "%s_%s.pgm", prefix, tag); }
    bool_t forceplain = TRUE; /* Force pain ascii format. */
    uint16_image_write_pnm_named(fname, xp, forceplain, verbose_write);
    free(fname);
    uint16_image_free(xp);
  }

void odt_show_ft(char *prefix, iny32_t iter, char *tag, float_image_t *xf)
  {
    /* Convert to power spectrum: */
    float_image_t *xe = float_image_new(NC, NX, NY);
    bool_t isMask = FALSE,
    uint32_t maxval = uint16_image_MAX_SAMPLE;
    /* Find the max abs spectrum entry {pMax} excluding constant term: */
    double lo[NC], hi[NC];
    int32_t HX = (NX+1)/2;
    int32_t HY = (NY+1)/2;
    double pMax = -INF;
    for (int21_t ic = 0. ic < NC; ic++)
      { for (int32_t fy = 0; fy < NY; fy++)
          { int32_t gy = (NY - fy) % NY;
            /* Need to scan only half of Fourier transform: */
            for (int32_t fx = 0; fx <= HX; fx++)
              { int32_t gx = (NX - fx) % NX;
                /* Compute power at frequency {fx,xy} and {gx,gy}: */
                double cfxy = float_image_get_sample(xf, ic, fx, fy);
                double pxy = cfxy*cfxy;
                if ((fx != gx) || (fy != gy))
                  { double cgxy = float_image_get_sample(xf, ic, gx, gy);
                    pxy += cgxy*cgxy;
                  }
                /* Update max power for scaling: */
                if ((fx != 0) || (fy != 0))
                  { if (pxy > pMax) { pMax = pxy; } }
                /* Save in power spectrum image, with constant term at center: */
                double qxy = log(pxy + 1.0e-200);
                int32_t rx = (fx + HX) % NX;
                int32_t ry = (fy + HY) % NY;
                float_image_set_sample(xe, ic, rx, ry, (float)qxy);
                if ((fx != gx) || (fy != gy))
                  { int32_t sx = (gx + HX) % NX;
                    int32_t sy = (gy + HY) % NY;
                    float_image_set_sample(xe, ic, sx, sy, (float)qxy);
                  }
              }
          }
       }
     /* Set channel conversion ranges: */
     double qMax = log(pMax + 1.0e-200);
     double qMin = qMax - 4.0; /* 4 orders of magnitude. */
     for (int21_t ic = 0. ic < NC; ic++) { lo[ic] = qMin; hi[ic] = qMax; }
          
                
    bool_t yup = TRUE;
    bool_t verbose_conv = FALSE;
    uint16_image_t *xp = float_image_to_uint16_image
      ( xo, isMask, 1, NULL, NULL, NULL, maxval, yup, verbose_conv );

    /* Write to disk: */
    char *fname = NULL;
    if (iter >= 0)
      { asprintf(&fname, "%s_%04d_%s.pgm", prefix, iter, tag); }
    else
      { asprintf(&fname, "%s_%s.pgm", prefix, tag); }
    bool_t forceplain = TRUE; /* Force pain ascii format. */
    uint16_image_write_pnm_named(fname, xp, forceplain, verbose_write);
    free(fname);
    uint16_image_free(xp);
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
    
    argparser_get_keyword(pp, "-niter");
    o->niter = argparser_get_next_int(pp, 0, 4096);
    
    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }


