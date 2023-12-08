#define PROG_NAME "geostereo"
#define PROG_DESC "extract stereo information from two anymaps"
#define PROG_VERS "1.0"

/* Copyright © 2003 by the State University of Campinas (UNICAMP).
** See the copyright, authorship, and warranty notice at end of file.
** Last edited on 2017-06-25 18:43:35 by stolfilocal
*/

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -mindisp {MINDISP} -maxdisp {MAXDISP} \\\n" \
  "  [ -window {WIDTH} {HEIGHT} ] \\\n" \
  "  [ -output {PREFIX} | [ -dispmap {DISPFILE} ] [ -scoremap {SCOREFILE} ] ] \\\n" \
  "  {PNMFILE1} {PNMFILE2}"
  
#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads two images of a 3D scene and computes" \
  " their stereo discrepancy tomogram {D}.  It also outputs the most likely" \
  " displacement map {dopt} and its confidence map {copt}.\n" \
  "\n" \
  "  In the discrepancy tomogram, each voxel is associated to" \
  " a triple {(x,y,d)} where {(x,y)} is" \
  " point of the image, and {d} is a displacement between" \
  " {MINDISP} and {MAXDISP}.  The voxel's value" \
  " is the discrepancy {D} between the two images in" \
  " the neighborhood of the points {(x+d,y)} and {(x-d,y)}.\n" \
  "\n" \
  "  In the {dopt} and {copt} maps, each pixel is associated" \
  " to a point {(x,y)} of the scene.  The most likely displacement" \
  " {dopt(x,y)} is the displacement {d} that minimizes" \
  " {D(x,y,d)}.  The  confidence {sopt[x,y]} is the value of {D(x,y,d)}" \
  " for that {d} !!! FINISH. \n" \
  "\n" \
  "DISCREPANCY FORMULA\n" \
  "  Presently, the discrepancy is the Hann-weighted root-mean-square" \
  " pixel difference within windows of specified size centered" \
  " at the two points.\n" \
  "\n" \
  "OUTPUT FORMAT\n" \
  "  !!! To be written.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -maxdisp {MAXDISP}\n" \
  "    Specifies the maximum relative displacement to consider, in pixels.  Required.\n" \
  "\n" \
  "  -mindisp {MINDISP} \n" \
  "    Specifies the minimum relative displacement to consider, in pixels.  Required.\n" \
  "\n" \
  "  -window {WIDTH} {HEIGHT}\n" \
  "    Specifies the width and height of the window used for" \
  " local matching.  Both numbers must be odd, and default to 3.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnm(5).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2003-02-07 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP).\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2003-02-07 created (J.Stolfi).\n" \
  "  2017-06-24 Changed to {float_image_t} throughout (J.Stolfi)."
  
/* TO DO: !!! Implement the tomogram output. */
#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <values.h>
#include <string.h>
#include <math.h>
 
#include <jsstring.h>
#include <argparser.h>

#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>
#include <float_image_geostereo.h>
#include <sample_conv.h>

???
void float_image_geostereo_normalize_samples
  ( float smp[], 
    uint32_t nwx, 
    uint32_t nwy, 
    double wt[],
    uint32_t NC
  );
  /* Independently normalizes each channel of the given samples {smp[]}
    to have mean 0 and unit variance. Ignores {NAN} samples. If all
    samples are equal, sets them all to 0.
    
    Assumes that the samples are stored in {smp} as
    described under {float_image_geostereo_get_samples}. */

???
void float_image_geostereo_normalize_samples
  ( float smp[], 
    uint32_t nwx, 
    double wtx[],
    uint32_t nwy, 
    double wty[],
    uint32_t NC
  )
  { int32_t i, c;
    for (c = 0; c < NC; c++)
      { /* Compute average {avg}: */
        double sum_ws = 0.0;
        double sum_w = 1.0e-200; /* To avoid division by zero if all are {NAN}. */
        /* Shift so that mean of valid pixels is 0: */
        s = 0.0; nok = 0;
        for (i = c; i < npix; i += NC) 
          { double wi = w[i]; if (! isnan(wi)) { s += wi; nok++; } }
        if (nok == 0) { return; }
        s /= (double)nok;
        for (i = c; i < npix; i += NC)
          { if (! isnan(w[i])) { w[i] = (float)(w[i] - s); } }
        /* Scale so that variance of valid pixels is 1: */
        s = 0.0;
        for (i = c; i < npix; i += NC)
          { double wi = w[i]; if (! isnan(wi)) { s += wi*wi; } }
        s = sqrt(s/(double)nok);
        if (s == 0.0) { return; }
        for (i = c; i < npix; i += NC) 
          { if (! isnan(w[i])) { w[i] = (float)(w[i]/s); } }
      }
  }




/* ====================================================================== */
/* INTERNAL PROTOTYPES */

typedef struct geostereo_options_t
  { char* inf1;      /* Filename of left eye view image. */
    char* inf2;      /* Filename of right eye view image. */
    double mindisp;  /* Minimum displacement to try (pixels). */
    double maxdisp;  /* Maximum displacement to try (pixels). */
    int32_t nscales; /* Number of scales in multiscale search (0 = single-scale). */
    int32_t wx;      /* Window width (must be odd). */
    int32_t wy;      /* Window height (must be odd). */
    char* otfd;      /* Filename for the displacement map image. */
    char* otfs;      /* Filename for the score map image. */
  } geostereo_options_t;

int32_t Log2(int32_t x); 
  /* {ceil(log_2(x))} */

int32_t main(int32_t argc, char **argv);

geostereo_options_t *geostereo_parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
  
void geostereo_read_images(geostereo_options_t *o, float_image_t **img1P, float_image_t **img2P);
  /* Reads the two imput images from files named {o->inf1} and {o->inf2}.
    Checks if they have the same size and channel count (1 for grayscale,
    3 for RGB).  Returns pointers to them in {*img1P} and {*img2P}. */

void geostereo_compute_displacement_map
  ( float_image_t *img1,  /* Image 1. */
    float_image_t *img2,  /* Image 2. */
    geostereo_options_t *o,         /* Command-line options. */
    float_image_t **imgdP, /* (OUT) Dispmap image. */
    float_image_t **imgsP  /* (OUT) Scoremap image. */
  );
  /* Stores in {*imgdP} a pointer to a newly allocated grayscale image
    {imgd} where {imgd[y,x]} is the best displacement found centered at
    pixel {x,y}.
    
    Also stores in {*imgsP} a ponter to a newly allocated image {imgs}
    such that {imgs[y,x]} the score (mismatch) corresponding to that
    diplacement. */

void geostereo_write_image(float_image_t *img, char *label, char *fname, bool_t isMask, double vMin, double vmax);
  /* Writes the image {img} to a PGM or PPM file file named {fname}.  
    Samples are linearly scaled and quantized from {[vmin _ vmax]} to {0..65535}.
    For the {isMask} parameter, see {sample_conv_quantize}. */

/* ====================================================================== */
/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { 
    /* Parse command line options: */
    geostereo_options_t *o = geostereo_parse_options(argc, argv);
  
    /* Read input images {img1,img2} and sizes {NC,NX,NY}: */
    float_image_t *img1, *img2;
    geostereo_read_images(o, &img1, &img2);
  
    /* Feature match: */
    float_image_t *imgd, *imgs;
    geostereo_compute_displacement_map(img1, img2, o, &imgd, &imgs);
      
    if (normalize)
      { /* Put NANs in both pl: */
        float_image_geostereo_copy_nans(smp1, smp2, nwx, nwy, NC);
        /* Normalize samples for zero mean and unit variance: */
        float_image_geostereo_normalize_samples(smp1, nwx, nwy, wt, NC);
        float_image_geostereo_normalize_samples(smp2, nwx, nwy, wt, NC);
        if (debug) 
          { fprintf(stderr, "    normalized samples:\n");
            float_image_geostereo_debug_window(smp1, nwx, nwy, NC); 
            float_image_geostereo_debug_window(smp2, nwx, nwy, NC);
          }
      }
      


    /* Write displacement map: */
    bool_t disp_isMask = FALSE; /* Smooth sample distr. */
    geostereo_write_image(imgd, "displacement", o->otfd, disp_isMask, o->mindisp, o->maxdisp);

    bool_t score_isMask = TRUE; /* Values 0 and 1 are special. */
    geostereo_write_image(imgs, "score", o->otfs, score_isMask, 0.0, 1.0);
    return 0;
  }
    



void geostereo_write_image(float_image_t *img, char *label, char *fname, bool_t isMask, double vMin, double vmax)
  { if ((fname == NULL) || ((*fname) == '\000')) { return; }
    fprintf(stderr, "writing the %s map...\n", label);
    bool_t yup = FALSE;
    bool_t warn_open = TRUE;
    bool_t debug_conv = FALSE;
    float_image_write_pnm_named(fname, img, isMask, 1.0, 0.0, yup, warn_open, debug_conv);
  }

void geostereo_read_images(geostereo_options_t *o, float_image_t **img1P, float_image_t **img2P) 
  { 
    bool_t isMask = FALSE; /* Smooth sample distr. */
    double gamma = sample_conv_BT709_DEC_GAMMA;
    double bias = sample_conv_BT709_BIAS;
    bool_t yup = FALSE;
    bool_t open_warn = TRUE;
    bool_t debug_conv = FALSE;
    float_image_t *img1 = float_image_read_pnm_named(o->inf1, isMask, gamma, bias, yup, open_warn, debug_conv);
    float_image_t *img2 = float_image_read_pnm_named(o->inf2, isMask, gamma, bias, yup, open_warn, debug_conv);
    
    /* Return images: */
    (*img1P) = img1;
    (*img2P) = img2;
  }

void geostereo_compute_displacement_map
  ( float_image_t *img1,    /* Image 1. */
    float_image_t *img2,    /* Image 2. */
    geostereo_options_t *o, /* Command-line options. */
    float_image_t **imgdP,  /* (OUT) Dispmap image. */
    float_image_t **imgsP   /* (OUT) Scoremap image. */
  )
  {
    /* Check image sizes: */
    int32_t NC, NX, NY;
    float_image_get_size(img1 ,&NC, &NX, &NY);
    float_image_check_size(img2, NC, NX, NY);

    /* Image pyramids: */
    
    /* Compute floating-point displacement map: */
    if (o->nscales > 0) 
      { fprintf(stderr, "using multiscale feature matching with %d levels\n", o->nscales);  } 
    else
      { fprintf(stderr, "using single-scale feature matching\n");  }
    float_image_t *imgd;  /* Displacement map. */
    float_image_t *imgs;  /* Score map. */
    int32_t ncands = 1; /* For now, keep only the best candidate. */
    float_image_geostereo_multiscale_displacement_map
      ( img1, img2, 
        /* nscales: */ o->nscales, 
        /* ncands: */ ncands,
        /* rx,ry: */ o->wx/2, o->wy/2,
        /* dmin,dmax: */ o->mindisp, o->maxdisp,
        &imgd, &imgs
      );
    
    /* Check dimensions of returned images: */
    float_image_check_size(imgd, ncands, NX, NY);
    float_image_check_size(imgs, ncands, NX, NY);
    
    /* Return images: */
    (*imgdP) = imgd;
    (*imgsP) = imgs;
  }
  
int32_t Log2(int32_t x)
  {
    int32_t e = 0;
    if (x < 0) { x = -x; }
    assert(x > 0);
    while (x > 1) { e++; x /= 2; }
    return e;
  }
  
#define geostereo_MAX_MAXDISP (8*1024)
#define geostereo_MAX_NSCALES 20
#define geostereo_MAX_WINDOW_SIZE 255

geostereo_options_t *geostereo_parse_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
     
    geostereo_options_t *o = (geostereo_options_t *)notnull(malloc(sizeof(geostereo_options_t)), "no mem");

    if (argparser_keyword_present(pp, "-window"))
      { o->wx = (int32_t)argparser_get_next_int(pp, 1, geostereo_MAX_WINDOW_SIZE); 
        if ((o->wx % 2) != 1)
          { argparser_error(pp, "window width must be odd"); } 
        o->wy = (int32_t)argparser_get_next_int(pp, 1, geostereo_MAX_WINDOW_SIZE);
        if ((o->wy % 2) != 1)
          { argparser_error(pp, "window height must be odd"); } 
      }
    else
      { o->wx = 3; o->wy = 3; }
      
    argparser_get_keyword(pp, "-maxdisp");
    o->maxdisp = argparser_get_next_double(pp, -geostereo_MAX_MAXDISP, +geostereo_MAX_MAXDISP);
      
    argparser_get_keyword(pp, "-mindisp");
    o->mindisp = argparser_get_next_double(pp, -geostereo_MAX_MAXDISP, o->maxdisp);
      
    if (argparser_keyword_present(pp, "-nscales"))
      { o->nscales = (int32_t)argparser_get_next_int(pp, 0, geostereo_MAX_NSCALES); }
    else 
      { /* Use enough scales to reduce images to to 1x1: */
        int32_t wmin = (o->wx < o->wy ? o->wx : o->wy);
        o->nscales = Log2(wmin);
      }
      
    if (argparser_keyword_present(pp, "-output"))
      { /* Output file names derived from a common prefix: */
        char *out_prefix = argparser_get_next(pp);
        o->otfd = txtcat(out_prefix, "-disp.pgm"); 
        o->otfs = txtcat(out_prefix, "-score.pgm"); 
      }
    else
      { /* Output file names specified separately: */
        if (argparser_keyword_present(pp, "-dispmap"))
          { o->otfd = argparser_get_next(pp); }
        else 
          { o->otfd = NULL; }

        if (argparser_keyword_present(pp, "-scoremap"))
          { o->otfs = argparser_get_next(pp); }
        else 
          { o->otfs = NULL; }
      }
    
    if ((o->otfd == NULL) && (o->otfs == NULL))
      { fprintf(stderr, "what, no output files?\n"); }
    
    /* Get input file names: */
    argparser_skip_parsed(pp);
    o->inf1 = argparser_get_next(pp);
    o->inf2 = argparser_get_next(pp);
    
    /* Check for extraneous arguments: */
    argparser_finish(pp);
    
    return o;
  }

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 2003 by the State University of Campinas (UNICAMP).
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
