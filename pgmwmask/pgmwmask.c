#define PROG_NAME "pgmwmask"
#define PROG_DESC "create a PGM weight mask for convolution, filtering, etc."
#define PROG_VERS "1.0"

/* Copyright © 2006 the State University of Campinas (UNICAMP).
** See the authorship, rights and warranty notices in the PROG_INFO below.
** Last edited on 2017-06-30 01:04:16 by stolfilocal
*/

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [ -window\\\n" \
  "      { square {W} | rect {WX} {WY} \\\n" \
  "      | round {W} | oval {WX} {WY} \\\n" \
  "      } \\\n" \
  "  ] \\\n" \
  "  [ -weights\\\n" \
  "      { uniform \\\n" \
  "      | gaussian {SX} {SY} \\\n" \
  "      | power {SX} {SY} {PWR} \\\n" \
  "      } \\\n" \
  "  ] \\\n" \
  "  [ -trimming { step | quadratic | biquadratic } ] \\\n" \
  "  [ -self {WT_SELF} ] \\\n" \
  "  [ -maxval {OMAXVAL} ] \\\n" \
  "  [ -output { text | pgm | PGM } ] \\\n" \
  "  [ -verbose ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program writes to stdout a weight mask" \
  " useful for window-based image filtering, such as convolution," \
  " median, local normalization, etc..\n" \
  "\n" \
  "   The mask's window (the pixels with nonzero weight) may be either" \
  " rectangular or round.  All samples outside the region are set to zero." \
  " The samples iside this region are computed from certain built-in" \
  " functions.\n" \
  "\n" \
  "   The sample values are conceptually real numbers" \
  " between 0 and 1. In the output file, they are encoded as" \
  " integers between 0 and {MAXVAL}, in such a way that the" \
  " output value {IVAL} represents the real number {IVAL/MAXVAL}.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -window {WIN_KIND} {WPARAMS}.. \n" \
  "    Specifies the window type and its parameters.  The valid" \
  " alternatives for {WIN_KIND} and {WPARAMS} are:\n" \
  "\n" \
  "      rect {WX} {WY}\n" \
  "        The window is a rectangle with {WX} columns and {WY}" \
  " rows." \
  "\n" \
  "      square {W}\n" \
  "        The window is a square spanning {W} rows" \
  " and columns.  Equivalent to \"rect {W} {W}\"." \
  "\n" \
  "      oval {WX} {WY}\n" \
  "        The window is a digital ellipse spanning {WX} columns and {WY}" \
  " rows." \
  "\n" \
  "      round {W}\n" \
  "        The window is a digital circle spanning {W} rows" \
  " and columns.  Equivalent to \"oval {W} {W}\"." \
  "\n" \
  "  -weights {WT_KIND} {WT_PARAMS}.. \n" \
  "    Specifies a weight function {w(X,Y)} that will define" \
  " the value of every pixel within the window, as a function of the pixel center" \
  " coordinates {(X,Y)}.  The valid alternatives" \
  " for {WT_KIND} and {WT_PARAMS} are:\n" \
  "\n" \
  "      uniform\n" \
  "        The function {w} is the unit constant, {w(X,Y) = 1} for all {X,Y}.\n" \
  "\n" \
  "      gaussian {SX} {SY}\n" \
  "        The function {w} is a two-dimensional Gaussian bell," \
  " {w(X,Y) = G((X-XC)/SX)*G((Y-YC)/SY)} where" \
  " {G(z) = exp(-z^2/2)}.  The" \
  " scaling factors {SX,SY} must be positive.\n" \
  "\n" \
  "      power {SX} {SY} {PWR}\n" \
  "        The function {w} is a two-dimensional inverse-power-like" \
  " distribution, {w(X,Y) = hypot((X-XC)/SX, (Y-YC)/SY, 1)^{-PWR}}.  The" \
  " scaling factors {SX,SY} and the exponent {PWR} must be positive.\n" \
  "\n" \
  "    In any case, the distribution specified by the \"-weights\" option" \
  " is truncated to the interior of the region specified by the" \
  " \"-window\" option.  The default is \"-weights uniform\".\n" \
  "\n" \
  "  -trimming { step | quadratic | biquadratic }\n" \
  "    Specifies how the weight function specified by the \"weights\" option" \
  " is to be clipped to the region specified by the \"-window\" option.  The" \
  " nominal value {w(X,Y)} of each pixel will be multiplied by a" \
  " trimming function {t(X,Y)}, which is 1 at the center of the mask" \
  " and 0 outside the window region.  Inside a round or oval window, the trimming" \
  " mask is a function of the relative position {R} of the pixel along" \
  " the ray between the mask center and the the window's edge.  Inside" \
  " a square or rectangular window, the" \
  " relevant parameters are relative abscissa {XR = (X-XC)/XC} and the" \
  " relative ordinate {YR} = {(Y-YC)/YC}.  The choices for {t} inside" \
  " the window are:\n" \
  "\n" \
  "      step\n" \
  "        The trimming function is 1 for all pixels, meaning that the" \
  " values {w(X,Y)} are used without modification.\n" \
  "\n" \
  "      quadratic\n" \
  "        The trimming function {t(X,Y)} is {Q(R)} if the window is" \
  " oval, and {Q(XR)*Q(YR)} if it is rectangular, where {Q(R) = 1-R^2}.  This" \
  " function has value 1 and slope 0 at the" \
  " center of the mask, and falls to 0 with non-zero slope at the window's" \
  " edge.\n" \
  "\n" \
  "      biquadratic\n" \
  "        The trimming function {t(X,Y)} is {K(R)} if the window is" \
  " oval, and {K(XR)*K(YR)} if it is rectangular; where {K(R) = (1-R^2)^2}.  This" \
  " function has value 1 and slope 0 at the" \
  " center of the mask, and falls to 0 with slope 0 at the window's" \
  " edge.\n" \
  "\n" \
  "    The default is \"-trimming biquadratic\".  In any case, for" \
  " pixels that straddle the window's boundary the trimming" \
  " function {t(X,Y)} is adjusted so that it makes a smooth" \
  " transition to 0.\n" \
  "\n" \
  "  -self {WC} \n" \
  "    This option forces the weight of the central pixel" \
  " to be {WC}, irrespective of the value defined by" \
  " the \"-weight\" option.  This option" \
  " can be used only if the width and height are both odd.\n" \
  "\n" \
  "  -maxval {OMAXVAL} \n" \
  "    Defines the nominal maximum sample value for the output" \
  " image.  Must be an integer between 1" \
  " and " stringify(PNM_MAX_MAXVAL) ".  If not specified," \
  " defaults to " stringify(PNM_MAX_MAXVAL) ".\n" \
  "\n" \
  "  -output { text | pgm | PGM }\n" \
  "    This optional parameter specifies the output" \
  " format for the mask.  With the \"text\" option, each sample" \
  " is written to {stdout} in a separate line," \
  " containing the sample" \
  " indices {IX,IY} and the quantized sample value {VAL}, in" \
  " decimal, separated by spaces.  Pixel rows are" \
  " separated by blank lines.  With the \"pgm\" or \"PGM\" options," \
  " the output is in the \"raw\" (binary) variant of the PGM image" \
  " file format.  The default is \"-output pgm\"\n" \
  "\n" \
  "  -verbose \n" \
  "    If this flag is present, several diagnostic messages and mask statistics" \
  " are written to {stderr}.  The default is to supress such" \
  " diagnostics.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  pgmkernel(1), pamgauss(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 2006-11-18 by J. Stolfi (IC-UNICAMP).\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2006 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <values.h>

#include <jspnm.h> 
#include <jsfile.h> 
#include <uint16_image.h>
#include <float_image.h>
#include <float_image_mask.h>
#include <float_image_to_uint16_image.h>
#include <affirm.h> 
#include <argparser.h> 

/* DATA TYPES */

#define INF INFINITY

typedef enum { WT_KIND_UNIF, WT_KIND_GAUSS, WT_KIND_POWER } wt_kind_t;

typedef struct options_t 
  { /* Parameters of the trimming function {t}: */
    bool_t window_oval;  /* Window shape: FALSE = rectangle, TRUE = oval. */
    int wx, wy;          /* Window width and height (0 if not given). */
    /* Parameters of the weight function {w}: */
    wt_kind_t wt_kind;   /* Type of weight distribution. */
    double sx, sy;       /* Scales for {x} and {y} axes ({NAN} if not given). */
    double pwr;          /* Exponent for inverse-power distrib ({NAN} if not given). */
    /* Parameters of the trimming function {t}: */
    int trimming;        /* 0 = step, 1 = quadratic, 2 = biquadratic. */
    /* Other parameters: */
    double self;         /* Weight of center pixel. */
    uint16_t maxval; /* Output maxval. */
    bool_t verbose;      /* TRUE prints diagnostics and statistics. */
    bool_t output_pgm;   /* FALSE = text format, TRUE = PGM format. */
  } options_t;
  /* Arguments parsed from the command line. */

/* INTERNAL PROTOTYPES */

int main(int argc, char* argv[]);
options_t *parse_options(int argc, char **argv);
float_image_t *create_weight_mask(options_t *o);
  /* Creates a weight image according to the specs in {o}. */
    
void print_weights(FILE *wr, uint16_image_t *imsk);
  /* Prints the weights to {wr} in the \"text\" format, etc. */

/* IMPLEMENTATIONS */

int main(int argc, char* argv[])
  { /* Command line arguments: */
    options_t *o = parse_options(argc, argv);

    float_image_t *fmsk = create_weight_mask(o);
    
    if (o->verbose) 
      { /* Statistics: */
        float_image_mask_stats_t S = float_image_mask_stats_get(fmsk, 0);
        fprintf(stderr, "\n");
        fprintf(stderr, "mask properties:\n");
        float_image_mask_stats_print(stderr, &S);
        fprintf(stderr, "\n");
      }
    
    /* Write the mask so that 0 and {o->maxval} represent half-width intervals. */
    bool_t isMask = TRUE;
    uint16_image_t *imsk = float_image_to_uint16_image(fmsk, isMask, 1, NULL, NULL, NULL, o->maxval, TRUE, o->verbose );
    if (o->output_pgm) 
      { bool_t forceplain = FALSE;
        uint16_image_write_pnm_file(stdout, imsk, forceplain, o->verbose); }
    else 
      { print_weights(stdout, imsk); }
    
    fflush(stdout);
    return 0;
  }

float_image_t *create_weight_mask(options_t *o)
  {
    /* Get the window dims and half-dims: */
    int wx = o->wx; assert(wx != 0);
    int wy = o->wy; assert(wy != 0);
    
    /* Create mask image: */
    float_image_t *fmsk = float_image_new(1, wx, wy);
    
    /* Fill it with the requested window trimming function: */
    float_image_mask_window(fmsk, 0, o->trimming, o->window_oval);
    
    /* Multiply the mask by the requested weight function: */
    switch(o->wt_kind)
      { case WT_KIND_UNIF:
          /* Just the window trimming function: */
          break;
        case WT_KIND_GAUSS:
          float_image_mask_mul_gauss(fmsk, 0, o->sx, o->sy);
          break;
        case WT_KIND_POWER:
          float_image_mask_mul_power(fmsk, 0, o->sx, o->sy, o->pwr);
          break;
        default:
          assert(FALSE); /* Invalid {o->wt_kind}. */
      }
    
    if (! isnan(o->self)) 
      { /* Set center pixel to specified value: */
        assert(wx % 2 == 1);
        assert(wy % 2 == 1);
        float_image_set_sample(fmsk, 0, wx/2, wy/2, o->self); 
      }
    
    return fmsk;
  }

void print_weights(FILE *wr, uint16_image_t *imsk)
  {
    /* Get the window dimensions: */
    assert (imsk->chns == 1); 
    int wx = imsk->cols;
    int wy = imsk->rows;
    
    fprintf(wr, "# columns = %d\n", wx);
    fprintf(wr, "# rows = %d\n", wy);

    /* Sample scaling factor: */
    double mv = (double)imsk->maxval;

    int x, y;
    for (y = 0; y < wy; y++)
      { for (x = 0; x < wx; x++)
          { uint16_t smp = uint16_image_get_sample(imsk, 0, x, y);
            fprintf(wr, "%4d %4d %5u %8.6f\n", x, y, smp, smp/mv);
          }
        fprintf(wr, "\n");
      }
  }

#define MAX_WINDOW_SIZE INT_MAX
  /* Maximum window height or width. */

options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
    
    /* Get keyword arguments: */
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = argparser_get_next_int(pp, 0, PNM_FILE_MAX_MAXVAL); }
    else
      { /* Use input maxval: */
        o->maxval = PNM_FILE_MAX_MAXVAL;
      }
    
    /* Window kind and parameters: */
    if (argparser_keyword_present(pp, "-window"))
      { int dims = 0; /* Number of dims required for this window kind. */
        if (argparser_keyword_present_next(pp, "rect"))
          { o->window_oval = FALSE; dims = 2; }
        else if (argparser_keyword_present_next(pp, "square"))
          { o->window_oval = FALSE; dims = 1; }
        else if (argparser_keyword_present_next(pp, "oval"))
          { o->window_oval = TRUE; dims = 2; }
        else if (argparser_keyword_present_next(pp, "round"))
          { o->window_oval = TRUE; dims = 1; }
        else
          { argparser_error(pp, "invalid window kind"); }
        /* Parse window dimensions, if appropriate: */
        o->wx = argparser_get_next_int(pp, 1, MAX_WINDOW_SIZE); 
        if (dims > 1)
          { o->wy = argparser_get_next_int(pp, 1, MAX_WINDOW_SIZE); }
        else
          { o->wy = o->wx; }
      }
    else
      { o->window_oval = TRUE;
        o->wx = 5; o->wy = 5;
      }
      
    if (argparser_keyword_present(pp, "-self"))
      { o->self = argparser_get_next_double(pp, 0.0, 1.0e+5);
        if ((o->wx % 2 != 1) || (o->wy % 2 != 1))
          { pnm_error("window width and height must be odd to use \"-self\""); }
      }
    else
      { o->self = NAN; }
    
    /* Weight distribution kind and parameters: */
    o->pwr = NAN;
    o->sx = NAN;
    o->sy = NAN;
    if (argparser_keyword_present(pp, "-weights"))
      { if (argparser_keyword_present_next(pp, "uniform"))
          { /* Uniform weights: */
            o->wt_kind = WT_KIND_UNIF;
          }
        else if (argparser_keyword_present_next(pp, "gaussian"))
          { /* Gaussian weight distribution: */
            o->wt_kind = WT_KIND_GAUSS;
            o->sx = argparser_get_next_double(pp, 1.0e-3, 1.0e+4);
            o->sy = argparser_get_next_double(pp, 1.0e-3, 1.0e+4);
          }
        else if (argparser_keyword_present_next(pp, "power"))
          { /* Inverse power law: */
            o->wt_kind = WT_KIND_POWER;
            o->sx = argparser_get_next_double(pp, 1.0, 1000.0);
            o->sy = argparser_get_next_double(pp, 1.0, 1000.0);
            o->pwr = argparser_get_next_double(pp, 0.1, 1000.0);
          }
        else
          { argparser_error(pp, "invalid weight function"); }
      }
    else
      { /* Default is uniform weights: */
        o->wt_kind = WT_KIND_UNIF;
      }
    
    /* Trimming function kind and parameters: */
    if (argparser_keyword_present(pp, "-trimming"))
      { if (argparser_keyword_present_next(pp, "step"))
          { o->trimming = 0; }
        else if (argparser_keyword_present_next(pp, "quadratic"))
          { o->trimming = 1; }
        else if (argparser_keyword_present_next(pp, "biquadratic"))
          { o->trimming = 2; }
        else
          { argparser_error(pp, "invalid trimming option"); }
      }
    else
      { /* Default is biquadratic trimming: */
        o->trimming = 2;
      }
    
    /* Show weights: */
    if (argparser_keyword_present(pp, "-output"))
      { if (argparser_keyword_present_next(pp, "text"))
          { o->output_pgm = FALSE; }
        else if (argparser_keyword_present_next(pp, "pgm"))
          { o->output_pgm = TRUE; }
        else if (argparser_keyword_present_next(pp, "PGM"))
          { o->output_pgm = TRUE; }
        else
          { argparser_error(pp, "invalid output_pgm option"); }
      }
    else
      { o->output_pgm = TRUE; }
    
    /* Diagnostics and statistics: */
    o->verbose = argparser_keyword_present(pp, "-verbose");
    
    /* Skip to positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }
