#define PROG_NAME "pnmfftfilter"
#define PROG_DESC "Fourier-based image filter"
#define PROG_VERS "1.0"

/* Copyright © 2008 by the State University of Campinas (UNICAMP).
** See the copyright, authorship, and warranty notice at end of file.
** Last edited on 2017-06-30 01:06:01 by stolfilocal
*/

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  { -pass | -kill } \\\n" \
  "  [ -from {W_MIN} [ {W_MIN_Y} ] ] \\\n" \
  "  [ -to {W_MAX} [ {W_MAX_Y} ] ] ] \\\n" \
  "  [ -range { {V_LO} {V_HI} | AUTO ] \\\n" \
  "  [ -maxval {MV_OUT} ] \\\n" \
  "  [ -verbose ] \\\n" \
  "  [ {PNMFILE_IN} ]"

#define PROG_INFO_DESC \
  "  The program reads the PNM image {PNMFILE_IN} and applies to" \
  " it a spatial band-pass or band-kill filter.\n" \
  "\n" \
  "  The band-pass filter produces the difference {F_UP - F_LO}" \
  " between two low-pass Gaussian filters {F_LO} and {F_UP}, with" \
  " characteristic spatial wavelengths {W_MAX} and {W_MIN}," \
  " respectively.  As special cases, the filter {F_UP} may be the" \
  " identity filter, and {F_LO} may be the zero filter." \
  "\n" \
  "  Thus, if {W[0],.. W[N]} is a list of {N+1} wavelengths, in" \
  " increasing order, with {W[0] = 0} and {W[N] = INF}, then" \
  " the images produced by \"" PROG_NAME " -from {W[i]} -to {W[i+1]}\"," \
  " with {i} in {0..N-1}, can be added together to reproduce the" \
  " original image (apart from the sample shift and scale introduced by" \
  " the conversion to the PBM format).\n" \
  "\n" \
  "  If the argument {PNMFILE_IN} is omitted or is \"-\", the" \
  " program reads the input image from {stdin}.  Each" \
  " sample {V_IN} is converted to a floating-point value in {[0_1]}" \
  " by {sample_conv_floatize(V_IN,MV_IN,FALSE,0.0,1.0,...)}." \
  "\n" \
  "  The output image is always written to {stdout}.  Each" \
  " sample {V_OUT} is rescaled and quantized" \
  " by {sample_conv_quantize(V_OUT,MV_OUT,FALSE,V_LO,V_HI,...)}."

#define PROG_INFO_OPTS \
  "  -pass\n" \
  "  -kill\n" \
  "    These options select between a band-pass filter (\"-pass\") or" \
  " a band-kill filter (\"-kill\").  If the \"-pass\" option attenuates" \
  " a Fourier component by some factor {W}, the \"-kill\" option attenuates" \
  " it by {1-W}.  At most one of the two options, \"-pass\" or \"-kill\"," \
  " must be specified.  The default is \"-pass\".\n" \
  "\n" \
  "  -from {W_MIN}\n" \
  "  -from {W_MIN_X} {W_MIN_Y}\n" \
  "    This option specifies the lower characteristic" \
  " wavelength (hence the upper characteristic frequency)" \
  " of the spatial filter.  In \"-pass\" mode" \
  " Fourier components with wavelengths below {W_MIN} pixels are" \
  " supressed.  If {W_MIN} is zero or negative, the" \
  " smallest-scale (highest-frequency) details are preserved, and" \
  " the program performs a high-frequency-pass instead of band-pass filter.  If" \
  " two values are given, they refer to the horizontal and vertical axes," \
  " respectively.  The wavelength(s) may be fractional.  The default" \
  " is \"-from 0\".\n" \
  "\n" \
  "  -to {W_MAX}\n" \
  "  -to {W_MAX_X} {W_MAX_Y}\n" \
  "    This option specifies the upper characteristic" \
  " wavelength (lower chracatteristic frequency) of the spatial" \
  " filter.  In \"-pass\" mode Fourier components with" \
  " wavelengths above {W_MAX} pixels are supressed.  If {W_MAX} is" \
  " finite but much larger than the image dimensions, the" \
  " program will remove only the zero" \
  " frequency (mean-value) term.  The value may also" \
  " be \"INF\", \"Inf\", \"inf\", \"+oo\" or \"oo\"," \
  " meaning infinity, in which case the program performs" \
  " a low-frequency-pass filter that preserve the mean value too.   If" \
  " two values are given, they refer to the horizontal and" \
  " vertical axes, respectively.  The wavelength(s) may be" \
  " fractional.   The default is \"-to +oo\".\n" \
  "\n" \
  "  -range {V_LO} {V_HI}\n" \
  "  -range AUTO\n" \
  "    This option specifies that all samples in the filtered" \
  " image must be rescaled before quantizing," \
  " so that the range {[V_LO _ V_HI]} is affinely mapped" \
  " to {[0 _ 1]}.   If the \"AUTO\" variant is used, or" \
  " if {V_LO} is equal to {V_HI}, sets {V_LO} and {V_HI} to" \
  " the maximum and minimum sample values in the filtered" \
  " image.  If this option is not given, assumes \"-range AUTO\".\n" \
  "\n" \
  "  -maxval {MV_OUT}\n" \
  "    Specifies {MV_OUT} as the maximum sample value for the" \
  " output image.  It must be an integer between 255 and 65535," \
  " inclusive. If not specified, it is set to 255 or to the" \
  " input image's {maxval}, whichever is larger.\n" \
  "\n" \
  "  -verbose\n" \
  "    If this option is present, the program prints out" \
  " global debugging information, such as input and output" \
  " image statistics." \
  
#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmscale(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created sep/2008 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  sep/2008 Created by adaptation of {pnmgtran.c}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2008 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <fftw3.h>

#include <jsfile.h>
#include <r2.h>
#include <uint16_image.h>
#include <float_image.h>
#include <float_image_hartley.h>
#include <float_image_filter.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <jspnm.h>
#include <argparser.h>

typedef struct options_t 
  { char *fname;         /* Input file name. */
    /* Characteristic wavelengths: */
    r2_t wMin;           /* Min wavelength in pixels. */
    r2_t wMax;           /* Max wavelength in pixels. */
    double v_lo;         /* Low output float value for scaling. */
    double v_hi;         /* High output float value for scaling. */
    bool_t complement;   /* TRUE for "-kill", FALSE for "-pass". */
    /* Output image attributes: */
    uint16_t maxval; /* Output maxval requested by user, or 0 if not given. */
    /* Debugging options: */
    bool_t verbose;      /* TRUE to print global statistics. */
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

options_t *get_options(int argc, char **argv);

r2_t parse_wavelength_argument(argparser_t *pp);
  /* Parses a wavelength argumetn, which may be either a float or a
    pair o floats. Either number may be "INF", "Inf", "inf" or "oo",
    possibly with a prefixed "+". If only one wavelength is given,
    sets both coordinates to that wavelength. */

float_image_t *read_image
  ( FILE *rd, 
    int *colsP, 
    int *rowsP, 
    int *chnsP, 
    uint16_t *maxvalP,
    bool_t verbose
  );
  /* Reads a PBM/PGM/PPM image file from {rd}, converts it to a float
    image with samples in the range [0_1]. Returns the relevant image
    data. If {verbose} is true, prints image statistics to
    {stderr}. */

void write_image
  ( FILE *wr, 
    float_image_t *fim, 
    double lo,
    double hi,
    uint16_t maxval,
    bool_t verbose
  );
  /* Writes the float image {fim} to {wr} as a PBM/PGM/PPM image file.
    Samples are converted from the range {[lo _ hi]} to {[0_1]}
    before being quantized to {0..maxval}.  However, if {lo == hi}, 
    sets those parameters to the min and max sample values instead.
    
    If {verbose} is true, prints image statistics to {stderr}. */
    
/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    /* Parse command line options: */
    options_t *o = get_options(argc, argv);

    /* Read input image, get dimensions: */
    int chns, cols, rows;
    uint16_t maxval_in;
    FILE *rd = open_read(o->fname, o->verbose);
    float_image_t *im_in = read_image(rd, &cols, &rows, &chns, &maxval_in, o->verbose);
    
    /* Compute Hartley transform of image: */
    float_image_t *im_ft = float_image_new(chns, cols, rows);
    float_image_hartley_transform(im_in, im_ft);
    
    /* Apply band filter: */
    float_image_filter_gaussian_band(im_ft, &(o->wMin), &(o->wMax), o->complement, o->verbose);
    
    /* Allocate output image: */
    float_image_t *im_ot = float_image_new(chns, cols, rows);
    
    /* Convert image back to original form: */
    float_image_hartley_transform(im_ft, im_ot);
    
    /* Choose output maxval: */
    uint16_t maxval_ot = (o->maxval > 0 ? o->maxval :(maxval_in < 255 ? 255 : maxval_in));

    /* Write the filtered image, scaled and quantized: */
    write_image(stdout, im_ot, o->v_lo, o->v_hi, maxval_ot, o->verbose);

    if (o->verbose) { fprintf(stderr, "done."); }
    return 0;
  }

float_image_t *read_image
  ( FILE *rd, 
    int *colsP, 
    int *rowsP, 
    int *chnsP, 
    uint16_t *maxvalP,
    bool_t verbose
  )
  { uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    bool_t isMask = FALSE; /* Assume smooth distr of sample values. */
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, TRUE, verbose);
    (*colsP) = pim->cols;
    (*rowsP) = pim->rows;
    (*chnsP) = pim->chns;
    (*maxvalP) = pim->maxval;
    uint16_image_free(pim);
    return fim;
  }

void write_image
  ( FILE *wr, 
    float_image_t *fim, 
    double lo,
    double hi,
    uint16_t maxval,
    bool_t verbose
  )
{ int chns = (int)fim->sz[0];
    double vLo[chns];
    double vHi[chns];
    int c;
    if (lo == hi)
      { float vMin = +INF, vMax = -INF;
        /* Find true range of samples over all channels: */
        for (c = 0; c < chns; c++) 
          { float_image_update_sample_range(fim, c, &vMin, &vMax); }
        for (c = 0; c < chns; c++) 
          { vLo[c] = vMin; vHi[c] = vMax; }
      }
    else
      { /* Use {[lo_hi]} for all channels: */
        for (c = 0; c < chns; c++) { vLo[c] = lo; vHi[c] = hi; }
      }
    bool_t isMask = FALSE; /* Assume smooth distr of sample values. */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, chns, vLo, vHi, NULL, maxval, TRUE, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    uint16_image_free(pim);
  }

#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

#define MAX_SIZE (32*1024)
  /* A limit on image size, to avoid humongous mallocs. */

options_t *get_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
     
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    if (argparser_keyword_present(pp, "-pass"))
      { o->complement = FALSE; }
    else if (argparser_keyword_present(pp, "-kill"))
      { o->complement = TRUE; }
    else
      { o->complement = FALSE; }
    
    if (argparser_keyword_present(pp, "-from"))
      { o->wMin = parse_wavelength_argument(pp); }
    else
      { o->wMin = (r2_t){{ 0, 0 }}; }
      
    if (argparser_keyword_present(pp, "-to"))
      { o->wMax = parse_wavelength_argument(pp); }
    else
      { o->wMax = (r2_t){{ +INF, +INF }}; }
      
    int ax;
    for (ax = 0; ax < 2; ax++)
      { if (o->wMin.c[ax] >= o->wMax.c[ax])
          { argparser_error(pp, "invalid wavelength band"); }
      }
   
    if (argparser_keyword_present(pp, "-range"))
      { if (argparser_keyword_present_next(pp, "AUTO"))
          { /* Automatic range: */
            o->v_lo = o->v_hi = 0;
          }
        else
          { /* User-specified range (auto if {v_lo >= v_hi}): */
            o->v_lo = argparser_get_next_double(pp, -INF, +INF);
            o->v_hi = argparser_get_next_double(pp, -INF, +INF);
          }
      }
    else
      { /* Automatic range: */
        o->v_lo = 0; o->v_hi = 0;
      }

    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (uint16_t)argparser_get_next_int(pp, 1, PNM_MAX_SAMPLE); }
    else
      { /* The default depends on the input image: */
        o->maxval = 0;
      }

    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Parse optional input file name: */
    argparser_skip_parsed(pp);
    if (argparser_next(pp) != NULL) 
      { o->fname = argparser_get_next(pp); }
    else
      { o->fname = "-"; }

    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }

r2_t parse_wavelength_argument(argparser_t *pp)
{
  r2_t w;
  int ax;
  for (ax = 0; ax < 2; ax++)
    { if 
        ( argparser_keyword_present_next(pp, "+INF") ||
          argparser_keyword_present_next(pp, "INF")  ||
          argparser_keyword_present_next(pp, "+Inf") ||
          argparser_keyword_present_next(pp, "Inf")  ||
          argparser_keyword_present_next(pp, "+inf") ||
          argparser_keyword_present_next(pp, "inf")  ||
          argparser_keyword_present_next(pp, "+oo")  ||
          argparser_keyword_present_next(pp, "oo")
        )
        { w.c[ax] = +INF; }
      else if (argparser_next_is_number(pp))
        { w.c[ax] = argparser_get_next_double(pp, 0.0, +INF); }
      else if (ax == 1)
        { w.c[ax] = w.c[0]; }
      else
        { argparser_error(pp, "expected a number or \"INF\""); }
    }
  return w;
 
}
