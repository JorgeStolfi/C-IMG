#define PROG_NAME "ppminksep"
#define PROG_DESC "ink analysis of digital images"
#define PROG_VERS "1.0"

/* Last edited on 2025-08-07 21:22:31 by stolfi */

/* Copyright © 2004 by the State University of Campinas (UNICAMP).*/
/* See the copyright, authorship, and warranty notice at end of file.*/

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    [ {-inGamma|-outGamma} {R_EXPO} {G_EXPO} {B_EXPO} ].. \\\n" \
  "    [ -logScale ] \\\n" \
  "    [ -layer {NAME} -color {INK_COLOR} [ -remap {NEW_COLOR} ] ]... \\\n" \
  "    [ -exceptions {NAME} ] \\\n" \
  "    [ -stretched {NAME} ] \\\n" \
  "    [ -shrink {SRINK_FRAC} ] \\\n" \
  "    [ -reveal ] \\\n" \
  "    [ -bgColor {BG_COLOR} ] \\\n" \
  "    [ -outPrefix {OUT_PREFIX} ] \\\n" \
  "    [ PNMFILE ]\n" \
  "  where each {NAME} is an output file identifier, \n" \
  "  each {*_EXPO} or {*_FRAC} is a value between 0 and 1, and\n" \
  "  each {*_COLOR} is " frgb_parse_color_HELP "."

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads a PPM (portable color pixmap) input file, and" \
  " analyzes it as a mixture of four given pigments (opaque" \
  " or transparent), in proportions that vary from pixel to" \
  " pixel.  Each pigment component is then written out as" \
  " a separate PPM file, overlaid on a background of arbitrary" \
  " color (which is specified by the \"-bgColor\" option).\n" \
  "\n" \
  "  If the \"-stretched\" option is given, the program will write out" \
  " also an image where the given colors are replaced by the corresponding" \
  " \"-map\" colors.  If these are not specified, the program" \
  " chooses suitable defaults.\n" \
  "\n" \
  "  The four pigment components can be interpreted as layers of partially" \
  " opaque paint stacked in a specific order.  The \"-reveal\" option" \
  " produces images that estimate the appearance of each layer if the" \
  " overaying layers were to be removed.\n" \
  "\n" \
  "  The program may also output an image showing the pixels" \
  " whose color cannot be explained as a mixture of the given" \
  " pigments. \n" \
  "\n" \
  PROG_INFO_COLOR_ANALYSIS "\n" \
  "\n" \
  PROG_INFO_REVEALING_COLORS "\n" \
  "\n" \
  PROG_INFO_TRANSPARENT_COLORS "\n" \
  "\n" \
  PROG_INFO_EXCEPTIONAL_PIXELS "\n" \
  "\n" \
  PROG_INFO_SHRINKING_COLORS "\n" \
  "\n" \
  PROG_INFO_GAMMA_CORRECTION "\n" \
  "\n" \
  "OPTIONS\n" \
  "  -layer {AMOUNT} {NAME} -color {INK_COLOR} [ -remap {NEW_COLOR} ] ]\n" \
  "    This option describes one ink layer, with nominal" \
  " color {INK_COLOR}.  The output file name will be \"{OUT_PREFIX}-{NAME}.ppm\".  If" \
  " the \"-remap\" option is present, the re-synthetized" \
  " image will use {NEW_COLOR} as the layer's ink" \
  " color; otherwise it will use {INK_COLOR}.  There" \
  " must be exactly four \"-layer\" arguments.\n" \
  "\n" \
  "    Each {*COLOR} argument consists of " frgb_parse_color_INFO "\n" \
  "\n" \
  "  The color components range beyween 0 and 1.\n" \
  "\n" \
  "  -outPrefix {OUT_PREFIX}\n" \
  "    This optional argument specifies a common prefix for all output" \
  " files. The file names will be \"{OUT_PREFIX}-{NAME}.{EXT}\" where" \
  " {NAME} is as specified in the \"-layer\", \"-exceptions\"," \
  " and \"-stretched\" options. If" \
  " omitted, it defaults to \"-outPrefix out/img\".\n" \
  "\n" \
  "  !!! FINISH THIS SECTION !!! \n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  \"Dr. Strangelove\" by Stanley Kubrick(1), if you haven't seen it already.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in 02/jul/2004 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2008-07-07 by J. Stolfi, IC-UNICAMP. Folded the manpage into the program.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2004 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_COLOR_ANALYSIS \
  "COLOR ANALYSIS\n" \
  "  Normally the program assumes that the components are" \
  " opaque pigments combined by mixing, so that the RGB value" \
  " of each pixel is a convex combination of the four given" \
  " RGB values.  The coefficients (which add to 1) represent" \
  " the percentage of the surface that is covered by each" \
  " pigment.  In the corresponding output images, these" \
  " same percentages of each color are overlaid on the" \
  " chosen output background color.  Note that if the" \
  " the paper or other backing material is partially" \
  " or totally visible in some parts of the input image, its" \
  " color must be specified as one of the four pigments."
  
#define PROG_INFO_REVEALING_COLORS \
  "REVEALING_COLORS\n" \
  "  If the \"-reveal\" option is given, the components" \
  " are interpreted as layers of opaque powdered pigments" \
  " stacked in the given order, with the first one at the" \
  " bottom (which is then interpreted as the backing" \
  " material).  The coefficient of a given layer in the" \
  " convex combination then represents the fraction of the" \
  " surface that is covered by its pigment and is not covered" \
  " by pigment from subsequent layers.  The program will" \
  " then try to recover the true coverage fraction of each layer," \
  " compensating for the effect of later layers.  That is," \
  " the adjusted coefficient assigned to layer {k} will then" \
  " represent the percentage of layers 1,2,... {k-1} that are" \
  " covered by the pigment of layer {k}.  Or, equivalently," \
  " the percentage of the surface that would be covered by" \
  " pigment {k} if layers {k+1}, {k+2}, ... 4 were not" \
  " present.  Thus, in particular, the first pigment" \
  " (backing material) will always get its coefficient" \
  " adjusted to coefficient 1.  Note that the adjusted" \
  " coefficient of a layer is indeterminate at a" \
  " pixel if subsequent layers cover that pixel completely."

#define PROG_INFO_TRANSPARENT_COLORS \
  "TRANSPARENT COLORS\n" \
  "  The convex combination model is appropriate for separating" \
  " layers of opaque pigments like tempera paints, or transparent" \
  " inks and dyes whose absorption spectra are disjoint.  However," \
  " transparent dyes with overlapping absorption spectra combine by" \
  " multiplication rather than convex combination, in which case the" \
  " convex combination analysis is not correct. For those cases," \
  " one should specify the \"-logScale\" option, which converts" \
  " all color data to logarithmic scale.  The multiplicative" \
  " behavior of the dyes then becomes additive, and the convex" \
  " analysis does apply.  The separated components are converted" \
  " back to linear scale before output, so the effect is like" \
  " chemically separating each dye and applying it on the" \
  " specified output background color."

#define PROG_INFO_EXCEPTIONAL_PIXELS \
  "EXCEPTIONAL PIXELS\n" \
  "  The given layer colors may be interpreted as corners" \
  " of a tetrahedron {T} in the RGB color space.  The colors" \
  " that lie inside {T} are those that can be obtained by mixing" \
  " positive amounts of the four given colors. Input pixels that" \
  " lie outside {T} are \"clipped\", i.e. replaced by the" \
  " ``nearest'' color in {T}, in a suitable sense, for the" \
  " purpose of the analysis. The \"-exceptions\" option writes" \
  " these ``exceptional'' pixels out to a fifth PPM image.  In this" \
  " image, pixels that were sucessfully analyzed as positive" \
  " combinations are painted with the output background color." \
  " specified output background color."

#define PROG_INFO_SHRINKING_COLORS \
  "SHRINKING_COLORS\n" \
  "  The \"-shrink\" option can be helpful when analyzing images" \
  " with many exceptional pixels.  It adjusts the computed" \
  " coefficients of the convex combination by making them" \
  " closer to 0.25.  The argument {Ns} of this option" \
  " controls the amount of this adjustment: 1.0 means no" \
  " change, 0.0 means make all coeffs equal to 0.25.  The" \
  " effect is equivalent to modifying the layer colors so" \
  " as to expand the tetrahedron {T} relative to its center" \
  " by a factor {1/Ns} (except that the output images will" \
  " still use the given colors). Thus, when {Ns} is less than 1, there" \
  " will be fewer exceptional pixels, but the components may not" \
  " be perfectly separated."

#define PROG_INFO_GAMMA_CORRECTION \
  "GAMMA CORRECTION\n" \
  "  By default, the program assumes that input and output" \
  " pixel values are proportional to light intensity.  However" \
  " images produced by most capture devices are encoded with a" \
  " power-like light transfer curve; in that case the exponent" \
  " of the curve (gamma) should be specified with the \"-inGamma\" option.  For" \
  " generality, the option takes three arguments -- one per channel -- although" \
  " they are usually equal. In that case one may also want to use" \
  " the \"-outGamma\" option, to have all output images encoded" \
  " in the same way."

#include <math.h>

#include <values.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>

#include <r4.h>
#include <r4x4.h>
#include <rmxn.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <frgb_inksep.h>
#include <jspnm.h>
#include <jsprintf.h>
#include <argparser.h>
#include <jsfile.h>
#include <uint16_image.h>
#include <sample_conv.h>
#include <sample_conv_gamma.h>

#define CHNS frgb_inksep_CHNS
#define LAYS frgb_inksep_LAYS

#define Black   (frgb_t){{0.0, 0.0, 0.0}}
#define White   (frgb_t){{1.0, 1.0, 1.0}}
#define Ones    (frgb_t){{1.0, 1.0, 1.0}}
#define NoColor (frgb_t){{NAN, NAN, NAN}}
#define LumDir  (frgb_t){{0.299f, 0.587f, 0.114f}}

#define MAX_PIXELS (16*256*256*256)   
#define MAX_SIZE MAX_PIXELS
  /* Maximum rows, columns, elements (to avoid absurd allocs). */

#define BT_ENC_BIAS sample_conv_gamma_BT709_BIAS
  /* Value of {bias} for {sample_conv_gamma} 
    that approximates the BT.709 encoding. */

#define MAXOUTFILES (LAYS+2)
  /* Output files are one per layer, plus one for exceptions and one for stretched: */

typedef struct layer_data_t
  { char *name;    /* Output file name. */
    frgb_t color;  /* Color of pigment on this layer. */
    frgb_t remap;  /* Color for re-inking. */
  } layer_data_t;

typedef struct options_t 
  { char *f_in_name;          /* Input filename, or "-" for stdin. */
    frgb_t inGamma_expo;           /* Expo of input image, per channel. */
    int32_t logScale;         /* TRUE to do computations in log scale. */
    layer_data_t layer[LAYS]; /* Specs for each color layer. */
    char *exceptions;         /* File name for exceptional pixels, or "". */
    char *stretched;          /* File name for stretched-color file, or "". */
    int32_t reveal;           /* TRUE to extrapolate hidden coverage factors. */
    double shrink;            /* Shrinking factor for barycentric coordinates. */
    frgb_t bgColor;           /* Background color to use for separated layers. */
    frgb_t outGamma_expo;          /* Expo of output image, per channel. */
    char *outPrefix;          /* Output filename prefix. */
  } options_t;
  
typedef struct color_scale_t
  { /* Gamma encoding exponents, per channel: */
    int32_t do_gamma;   /* TRUE if gamma correction is not trivial. */
    frgb_t expo;        /* Rel intensity is {(pixval/maxval)**expo}. */
    /* Parameters for input/output logscale conversion: */
    double eps;         /* Use linear scale for abs values below this. */
    double lmag;        /* Parameter for logscale. */
  } color_scale_t;   

#define TwoPi (2.0 * M_PI)

int32_t main_debug = TRUE;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);

FILE *pis_open_write_ppm(char *outPrefix, char *imgName, bool_t verbose);
  /* Opens the file "{outPrefix}-{imgName}.ppm" for writing.  
    If {verbose} is true, prints a message to {stderr}. */

options_t *get_options(int32_t argc, char **argv);

double parse_double(int32_t *argn, int32_t argc, char **argv, double lo, double hi);
  /* Parses the next command line argument {argv[*argn]} as a 
    floating-point number.  Complains if the number is not in {[lo _ hi]}.
    Increments {*argn}. */
    
int32_t parse_int(int32_t *argn, int32_t argc, char **argv, int32_t lo, int32_t hi);
  /* Parses the next command line argument {argv[*argn]} as an 
    integer number.  Complains if the number is not in {[lo .. hi]}.
    Increments {*argn}. */

int32_t eq_color(frgb_t a, frgb_t b);
  /* TRUE iff {a.c[i]} equals {b.c[i]} for {i} in {0..CHNS-1}. */
  
frgb_t max_color(frgb_t *ctr, frgb_t *dir);
  /* Returns the color {r = ctr + t*dir} where {t} is the largest
    positive real such that {r} is inside the unit cube. The vector
    {dir} should not be zero, and {ctr} must lie inside the unit
    cube. */

void print_color(FILE *f, char *pref, frgb_t *fv, char *suff);
  /* Print the given color to file {f}, surrounded by the given {pref} and
    {suff} strings. */
    
void fix_remap_colors(layer_data_t *layer, double y_min, double y_med, double y_max);
  /* Provides suitable defaults for {layer[k].remap}, if it is {NoColor}.
    The darkest and lightest colors will be mapped to greys of brightness
    {y_min} and {y_max}, repectively.  The other two will be mapped to 
    maximally saturated colors of brightness {y_med}. */
  
double brightness(frgb_t *c);
  /* Returns the perceptual brightness of {c}. */

/* COLOR SCALE CORRECTION */

void lut_in(frgb_t *p, color_scale_t *cs);
  /* Applies color scale conversions to input pixel value {p}. */

void lut_ot(frgb_t *p, color_scale_t *cs);
  /* Applies inverse color scale conversions to pixel value {p} for output. */

double log_in(double y, double eps, double lmag);
  /* Changes input channel data {y} to log scale, by the map
    {lmag*eps*(log(y/eps) + 1)} when {y>eps}, {lmag*y} when {y < eps}. 
    For negative {y}, uses absolute value and negates result. */

double log_ot(double y, double eps, double lmag);
  /* Undoes the effect of {log_in} on the logscale intensity {y}. */

color_scale_t make_color_scale(frgb_t *expo, int32_t logScale, double eps);
  /* Creates a {color_scale_t} record for the given gamma exponent triplet
    log scale option {logScale}, and logscale threshold {eps}. */

/* OPTION PARSING */

frgb_t parse_triplet(argparser_t *pp, double lo, double hi);
 /* Parses three consecutive numbers from the command line,
   checks whether each component lies in {[lo _ hi]}. */
  
layer_data_t parse_layer_data(argparser_t *pp);
  /* Parses the next few command line arguments, starting at
    {argv[*argn]}, as a paint layer description. Increments {*argn}. */ 
 
layer_data_t *alloc_layer_data_record(int32_t np);
  /* Allocates a field specs record and its internal tables,
    with enough space for {np} samples. */
    
void process_pixels
  ( options_t *o, 
    uint32_t rows, 
    uint32_t cols, 
    FILE *f_in, 
    uint16_t maxval_in,
    bool_t raw_in,
    bool_t bits_in,
    FILE **f_ot, 
    uint16_t maxval_ot, 
    bool_t raw_ot, 
    bool_t bits_ot
  );
  /* Reads body of input image from file {f}, computes output pixels,
    writes them to stdout. */

void print_int_pixel(FILE *f, char *pref, int32_t *iv, char *suff);
  /* Print the given pixel to file {f}, surrounded by the given {pref} and
    {suff} strings. */

/* DEBUGGING AIDS */

void debug_int_pixel(char *label, int32_t col, int32_t row, int32_t *iv, char *tail);
  /* Print the given pixel if the flag {main_debug} is true. */

void debug_color(char *label, int32_t col, int32_t row, frgb_t *fv, char *tail);
  /* Print the given color if the flag {main_debug} is true. */

/* IMPLEMENTATIONS */

#ifndef PGM_OVERALLMAXVAL
#define PGM_OVERALLMAXVAL PGM_MAXMAXVAL
#endif

#define INF INFINITY   

#define arg_error(N,A,MSG) \
  pnm_error("argument %d = %s: %s", (N), (A), (MSG))

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    /* Open input file: */
    FILE *f_in = open_read(o->f_in_name, TRUE);       /* Input file. */
    
    /* Format specs from input file: */
    uint32_t rows, cols, chns;
    uint16_t imaxval;
    bool_t iraw, ibits;
    pnm_format_t iformat;
    pnm_read_header(f_in, &cols, &rows, &chns, &imaxval, &iraw, &ibits, &iformat);
    if (chns != CHNS) { pnm_error("bad input file format - must be a color image"); }
      
    /* Output files are one per layer, plus one for exceptions and one for stretched: */
    int32_t otfs  = LAYS + 2;
    int32_t ixf_exceptions = LAYS;
    int32_t ixf_stretched = LAYS + 1;

    /* Format specs for output file: */
    uint16_t omaxval = (imaxval < 255 ? 255 : imaxval);
    bool_t oraw, obits;
    pnm_format_t oformat;
    pnm_choose_output_format(omaxval, chns, (! iraw), &oformat, &oraw, &obits);
    
    /* Open output files and write headers: */
    FILE *f_ot[otfs];  /* Output files: sep layers, exceptions, stretched. */
    for (int32_t k = 0;  k < LAYS; k++)
      { f_ot[k] = pis_open_write_ppm(o->outPrefix, o->layer[k].name, TRUE);
        pnm_write_header(f_ot[k], cols, rows, omaxval, oformat);
      }

    /* Open output file for exceptional pixels: */
    if (o->exceptions != NULL)
      { f_ot[ixf_exceptions] = pis_open_write_ppm(o->outPrefix, o->exceptions, TRUE);
        pnm_write_header(f_ot[ixf_exceptions], cols, rows, omaxval, oformat);
      }
    else
      { f_ot[ixf_exceptions] = NULL; }
      
    /* Open output file for stretched-color image: */
    if (o->stretched != NULL)
      { f_ot[ixf_stretched] = pis_open_write_ppm(o->outPrefix, o->stretched, TRUE);
        pnm_write_header(f_ot[ixf_stretched], cols, rows, omaxval, oformat);
      }
    else
      { f_ot[ixf_stretched] = NULL; }

    /* Process image: */
    
    process_pixels
      ( o, rows, cols,
        f_in, imaxval, iraw, ibits,
        f_ot, omaxval, oraw, obits
      );

    for (int32_t k = 0;  k < otfs; k++) { if (f_ot[k] != NULL) { fclose(f_ot[k]); } }
    if (f_in != stdin) { fclose(f_in); }
    pnm_message("done.");
    return 0;
  }
  
FILE *pis_open_write_ppm(char *outPrefix, char *imgName, bool_t verbose)
  { char *fname = jsprintf("%s-%s.ppm", outPrefix, imgName);
    FILE *wr = open_write(fname, verbose);
    free(fname);
    return wr;
  }

void process_pixels
  ( options_t *o, 
    uint32_t rows, 
    uint32_t cols, 
    FILE *f_in, 
    uint16_t maxval_in,
    bool_t raw_in,
    bool_t bits_in,
    FILE **f_ot, 
    uint16_t maxval_ot, 
    bool_t raw_ot, 
    bool_t bits_ot
  )
  {
    int32_t otfs = LAYS + 2; /* Number of output files. */
    int32_t ixf_exceptions = LAYS;     /* Index of exceptions layer. */
    int32_t ixf_stretched = LAYS + 1;  /* Index of stretched-color image. */
    
    /* Input and output scaling parameters: */
    double zero_in = 0, unit_in = (double)maxval_in;
    double scale_in = unit_in - zero_in;
    
    double zero_ot = 0, unit_ot = (double)maxval_ot;
    double scale_ot = unit_ot - zero_ot;
    
    /* Minimum step for floated input colors: */
    double eps = 1.0/scale_in;
    
    /* Input/output conversion parameters: */
    color_scale_t cs_in = make_color_scale(&(o->inGamma_expo), o->logScale, eps);
    color_scale_t cs_ot = make_color_scale(&(o->outGamma_expo), o->logScale, eps);

    /* Convert color options to internal scale: */
    for (int32_t i = 0;  i < LAYS; i++)  
      { if (main_debug) { fprintf(stderr, "layer %s colors:\n", o->layer[i].name); }
        lut_in(&(o->layer[i].color), &cs_in);
        lut_in(&(o->layer[i].remap), &cs_in);
      }
    if (main_debug) { fprintf(stderr, "bg color:\n"); }
    lut_in(&(o->bgColor), &cs_in);
    
    /* Provide suitable defaults for stretched image colors: */
    if (o->logScale)
      { fix_remap_colors(o->layer, 0.25, 0.50, 1.00); }
    else
      { fix_remap_colors(o->layer, 0.00, 0.55, 1.00); }

    /* Synthesis matrices (each row is an ink color, homogenized): */
    assert(LAYS == 4);
    r4x4_t mix_org; /* Syntesis matrix for original inks. */
    r4x4_t mix_new; /* Syntesis matrix for new inks. */
    for (int32_t i = 0;  i < 4; i++)
      { debug_color("lc", -1, -1, &(o->layer[i].color), "\n");
        debug_color("lr", -1, -1, &(o->layer[i].remap), "\n");
        for (int32_t j = 0;  j < 4; j++)
          { mix_org.c[i][j] = (j < CHNS ? o->layer[i].color.c[j] : 1.0);
            mix_new.c[i][j] = (j < CHNS ? o->layer[i].remap.c[j] : 1.0);
          }
      }
    
    if (main_debug)
      { rmxn_gen_print2
          ( stderr, 4, 4, (double*)&(mix_org.c), 4, (double*)&(mix_new.c), "%8.5f", 
            "mix_org,mix_new = [\n  ", "\n  ", "\n]\n", "[ ", " ", " ]", "  "
          );
      }
    
    /* Shrink the layer color simplex if so requested: */
    frgb_inksep_shrink_color_set(&mix_org,  o->shrink);
    
    /* Decomposition matrix: */
    r4x4_t sep;
    r4x4_inv(&mix_org, &sep);
    
    if (main_debug)
      { rmxn_gen_print3
          ( stderr, 4, 4, (double*)&(mix_org.c), 4, (double*)&(mix_new.c), 4, (double*)&(sep.c), "%8.5f", 
            "mix_org,mix_new,sep = [\n  ", "\n  ", "\n]\n", "[ ", " ", " ]", "  "
          );
      }
    
    /* Scanline buffer for input pixels. */
    uint16_t *smp_in = uint16_image_alloc_pixel_row(cols, CHNS);

    /* Scanline buffers for each output file. */
    uint16_t *smp_ot[otfs];
    
    for (int32_t k = 0;  k < otfs; k++)
      { smp_ot[k] = uint16_image_alloc_pixel_row(cols, CHNS); }
      
    /* Loop on pixels: */
    uint16_t *xp_in, *xp_ot[otfs];
    for (int32_t row = 0; row < rows; ++row)
      { pnm_read_pixels(f_in, smp_in, cols, CHNS, maxval_in, raw_in, bits_in);
        xp_in = &(smp_in[0]);
        for (int32_t k = 0;  k < otfs; k++) { xp_ot[k] = &(smp_ot[k][0]); }
        for(int32_t col = 0; col < cols; ++col)
          { main_debug = (col == 15) && (row == 15);
            /* Process input image. Assumes PPM_FORMAT or RPPM_FORMAT: */
            frgb_t fv;
            for (int32_t c = 0;  c < CHNS; c++)
              { fv.c[c] = (float)frgb_floatize(*xp_in, maxval_in, zero_in, scale_in);
                xp_in++;
              }
                
            debug_color("qi", col, row,&fv, "\n");

            /* Apply input gamma and/or log conversion: */
            lut_in(&fv, &cs_in);
            debug_color("ci", col, row,&fv, "\n");
           
            /* Compute the ink mixing coefficients {m[0..LAYS-1]}: */
            r4_t m;
            frgb_inksep_separate_layers(fv, &sep, &m);
            if (main_debug)
              { r4_gen_print(stderr, &m, "%8.5f", "m = [ ", " ", " ]\n"); }
            
            /* Compute ink-substituted color {sv}, clip to unit cube: */
            frgb_t sv = frgb_inksep_remap_color(&m, &mix_new);
            frgb_clip_rgb_towards_grey(&sv);
            debug_color("cs", col, row,&sv, "\n");
            
            /* Clip mixing coefficients to canonical simplex: */
            int32_t nbad; /* Number of negative mixing coefficients. */
            frgb_inksep_clip_to_simplex(&m, &nbad);
             
            if (o->reveal) 
              { /* Extrapolate actual coverage of lower layers: */
                frgb_inksep_reveal_colors(&m);
              }
           
            /* Compute color separations {pv[0..NLAYS-1]}: */
            frgb_t pv[LAYS];
            frgb_inksep_compute_separations(o->bgColor, &m, &mix_org, pv);
        
            /* Compute exceptional pixel color {ev}: */
            frgb_t ev = (nbad > 0 ? fv : o->bgColor);
            debug_color("ce", col, row,&sv, "\n");
            
            /* Generate output pixels: */
            for (int32_t k = 0;  k < otfs; k++)
              { 
                frgb_t gv;
                int32_t iv[CHNS];
                
                if (f_ot[k] != NULL)
                  { 
                    if (k < LAYS)
                      { /* Output layer {k} of color separation: */
                        gv = pv[k];
                      }
                    else if (k == ixf_exceptions)
                      { /* Output exceptional pixel image: */
                        gv = ev; 
                      }
                    else if (k == ixf_stretched)
                      { /* Output reinked image: */
                        gv = sv;
                      }
                    else
                      { pnm_error("bad file index"); }
                    debug_color("cg", col, row,&gv, "\n");

                    /* Apply output un-log and gamma conversions: */
                    lut_ot(&gv, &cs_ot);  
                    debug_color("co", col, row,&gv, "\n");

                    for (int32_t i = 0;  i < CHNS; i++) 
                      { iv[i] = frgb_quantize(gv.c[i], zero_ot, scale_ot, maxval_ot); }
                    debug_int_pixel("qo", col, row,iv, "\n");
                    
                    for (int32_t i = 0;  i < CHNS; i++) { *(xp_ot[k]) = (uint16_t)iv[i]; ++xp_ot[k]; }
                    if (main_debug) { fprintf(stderr, "\n\n"); }
                  }
              }
          }
        for (int32_t k = 0;  k < otfs; k++)
          { if (f_ot[k] != NULL)
              { pnm_write_pixels(f_ot[k], smp_ot[k], cols, CHNS, maxval_ot, raw_ot, bits_ot); }
          }
      }
  }
  
color_scale_t make_color_scale(frgb_t *expo, int32_t logScale, double eps)
  {
    color_scale_t cs;
    cs.expo = (*expo);
    cs.do_gamma = (! eq_color(*expo, Ones));
    if (logScale)
      { cs.eps = eps;
        cs.lmag = (eps < 1.0 ? 1.0/eps/(1 - log(eps)) : 1.0);
      }
    else
      { cs.eps = INF;
        cs.lmag = 1.0;
      }
    return cs;
  }

void lut_in(frgb_t *p, color_scale_t *cs)
  {
    debug_color("  lut_in p in", -1, -1, p, "\n");
    if (cs->do_gamma) 
      { for (int32_t i = 0;  i < CHNS; i++)
          { p->c[i] = sample_conv_gamma(p->c[i], cs->expo.c[i], BT_ENC_BIAS); }
      }
    if (cs->eps < INF)
      { for (int32_t i = 0;  i < CHNS; i++)
          { p->c[i] = (float)log_in(p->c[i], cs->eps, cs->lmag); }
      }
    debug_color("  lut_in p ot", -1, -1, p, "\n");
  }

void lut_ot(frgb_t *p, color_scale_t *cs)
  {
    if (cs->eps < INF)
      { for (int32_t i = 0;  i < CHNS; i++)
          { p->c[i] = (float)log_ot(p->c[i], cs->eps, cs->lmag); }
      }
    if (cs->do_gamma) 
      { for (int32_t i = 0;  i < CHNS; i++)
          { p->c[i] = sample_conv_gamma(p->c[i], 1/cs->expo.c[i], BT_ENC_BIAS); }
      }
  }

double log_in(double y, double eps, double lmag)
  {
    if ((eps < INF) && (! isnan(y)) && (y != 0.0) && (y != 1.0) && (y < INF))
      { double yabs = fabs(y);
        if (yabs > eps)
          { double ylog = eps*(log(yabs/eps) + 1); 
            y = (y < 0.0 ? -ylog : ylog);
          }
        y *= lmag;
      }
    return y;
  }

double log_ot(double y, double eps, double lmag)
  {
    if ((eps < INF) && (! isnan(y)) && (y != 0.0) && (y != 1.0) && (y < INF))
      { y /= lmag;
        double ylog = fabs(y);
        if (ylog > eps)
          { double yabs = eps*exp(ylog/eps - 1);
            y = (y < 0.0 ? -yabs : yabs);
          }
      }
    return y;
  }
      
int32_t eq_color(frgb_t a, frgb_t b)
  { for (int32_t ch = 0;  ch < CHNS; ch++)
      { float ac = a.c[ch];
        float bc = b.c[ch];
        if ((ac != bc) || (isnan(ac) != isnan(bc))) { return FALSE; }
      }
    return TRUE;
  }
  
double brightness(frgb_t *c)
  { return LumDir.c[0]*c->c[0] + LumDir.c[1]*c->c[1] + LumDir.c[2]*c->c[2];
  }

void fix_remap_colors(layer_data_t *layer, double y_min, double y_med, double y_max)
  {
    /* Assumes {CHNS = 3}. */
    if (CHNS != 3) { pnm_error("oops, CHNS not 3!"); }
    /* Check if any map-color needs to be fixed: */
    int32_t undef = 0;
    for (int32_t k = 0;  k < LAYS; k++)
      { undef += eq_color(layer[k].remap, NoColor); }
    if (undef == 0) { return; } 
    if (undef < LAYS)
      { pnm_error("option \"-remap\" must be given for all layers or for none"); }
    /* Find the lightest layer {kwh}, map to white: */
    int32_t kwh = -1; double ymax = -INF;
    for (int32_t k = 0;  k < LAYS; k++)
      { double y = brightness(&(layer[k].color));
        if (y >= ymax) { ymax = y; kwh = k; }
      }
    layer[kwh].remap = (frgb_t){{ (float)y_max, (float)y_max, (float)y_max }};
    /* Find the darkest layer {kbk} other than {kwh}, map to black: */
    int32_t kbk = -1; double ymin = INF;
    for (int32_t k = 0;  k < LAYS; k++)
      { double y = brightness(&(layer[k].color));
        if (y <= ymin) { ymin = y; kbk = k; }
      }
    layer[kbk].remap = (frgb_t){{ (float)y_min, (float)y_min, (float)y_min }};

    /* Find mean color: */
    frgb_t ctr = Black;
    for (int32_t k = 0;  k < LAYS; k++)
      { ctr.c[0] += layer[k].color.c[0];
        ctr.c[1] += layer[k].color.c[1];
        ctr.c[2] += layer[k].color.c[2];
      }
    ctr.c[0] /= LAYS;
    ctr.c[1] /= LAYS;
    ctr.c[2] /= LAYS;

    /* Find the other two layers {ka,kb}: */
    int32_t ka, kb;
    ka = 0; 
    if ((ka == kwh) || (ka == kbk)) { ka = (ka + 1) % LAYS; }
    if ((ka == kwh) || (ka == kbk)) { ka = (ka + 1) % LAYS; }
    kb = (ka + 1) % LAYS; 
    if ((kb == kwh) || (kb == kbk)) { kb = (kb + 1) % LAYS; }
    if ((kb == kwh) || (kb == kbk)) { kb = (kb + 1) % LAYS; }
    
    /* Find mean {m} those two colors, subtract from mean color: */
    frgb_t m;
    m.c[0] = (layer[ka].color.c[0] + layer[kb].color.c[0])/2 - ctr.c[0];
    m.c[1] = (layer[ka].color.c[1] + layer[kb].color.c[1])/2 - ctr.c[1];
    m.c[2] = (layer[ka].color.c[2] + layer[kb].color.c[2])/2 - ctr.c[2];
    
    /* Turn {m} into a pure chroma vector: */
    double ym = brightness(&m);
    m = frgb_shift(-ym, &m);
    /* Make sure {m} is not too small: */
    if (fabs(m.c[0]) + fabs(m.c[1]) + fabs(m.c[2]) < 0.001)
      { m = (frgb_t){{ LumDir.c[1], -LumDir.c[0], 0.0 }}; }
    /* Normalize {m} to unit length: */
    double sm = sqrt(m.c[0]*m.c[0] + m.c[1]*m.c[1] + m.c[2]*m.c[2]) + 1.0e-200;
    m = frgb_scale(1/sm, &m);
    /* Find a unit-length pure chroma vector {n} orthogonal to {m}: */
    frgb_t n = (frgb_t){{ 
      LumDir.c[2]*m.c[1] - LumDir.c[1]*m.c[2], 
      LumDir.c[0]*m.c[2] - LumDir.c[2]*m.c[0], 
      LumDir.c[1]*m.c[0] - LumDir.c[0]*m.c[1]
    }};
    double sn = sqrt(n.c[0]*n.c[0] + n.c[1]*n.c[1] + n.c[2]*n.c[2]) + 1.0e-200;
    n = frgb_scale(1/sn, &n);
    /* Compute two orthogonal pure chroma vector {va,vb} centered at {m}: */
    frgb_t va, vb;
    va.c[0] = m.c[0] + n.c[0]; vb.c[0] = m.c[0] - n.c[0];
    va.c[1] = m.c[1] + n.c[1]; vb.c[1] = m.c[1] - n.c[1];
    va.c[2] = m.c[2] + n.c[2]; vb.c[2] = m.c[2] - n.c[2];
    /* Map layers {ka,kb} to the max-sat colors of lum=0.5 and those chromas: */
    frgb_t gray = (frgb_t){{ (float)y_med, (float)y_med, (float)y_med }};
    layer[ka].remap = max_color(&gray, &va);
    layer[kb].remap = max_color(&gray, &vb);

    /* Report choices: */
    pnm_message("default colors for stretched image:");
    for (int32_t k = 0;  k < LAYS; k++)
     { fprintf(stderr, "  layer %d", k);
       print_color(stderr, " = ", &(layer[k].remap), "");
       fprintf(stderr, "\n");
     }
  }

frgb_t max_color(frgb_t *ctr, frgb_t *dir)
  {
    double t = +INF;
    for (int32_t i = 0;  i < CHNS; i++)
      { double ci = ctr->c[i], di = dir->c[i];
        double ti = (di == 0 ? +INF : (di > 0 ? (1.0 - ci)/di : (0.0 - ci)/di));
        if (ti < t) { t = ti; }
      }
    if (t == +INF) 
      { return *ctr; }
    else
      { frgb_t res;
        for (int32_t i = 0;  i < CHNS; i++)
          { res.c[i] = (float)(ctr->c[i] + t*dir->c[i]); }
        return res;
      }
  }

void debug_int_pixel(char *label, int32_t col, int32_t row, int32_t *iv, char *tail)
  { 
    if (main_debug) 
      { fprintf(stderr, "%s[%3d][%3d] = ", label, row, col);
        print_int_pixel(stderr, "( ", iv, " )");
        fprintf(stderr, "%s", tail);
      }
  }

void print_color(FILE *f, char *pref, frgb_t *fv, char *suff)
  { 
    fprintf(f, "%s", pref);
    for (int32_t k = 0;  k < CHNS; k++) 
      { fprintf(f, "%s%7.4f", (k == 0 ? "" : " "), fv->c[k]); }
    fprintf(f, "%s", suff);
  }

void print_int_pixel(FILE *f, char *pref, int32_t *iv, char *suff)
  { 
    fprintf(f, "%s", pref);
    for (int32_t k = 0;  k < CHNS; k++) 
      { fprintf(stderr, "%s%3d", (k == 0 ? "" : " "), iv[k]); }
    fprintf(f, "%s", suff);
  }

options_t *get_options(int32_t argc, char **argv)
  {
    options_t *o = (options_t *)malloc(sizeof(options_t));
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Parse keyword arguments: */
    if (argparser_keyword_present(pp, "-inGamma"))
      { o->inGamma_expo = parse_triplet(pp, 0.1, 10.0); } 
    else
      { o->inGamma_expo = (frgb_t){{1.0, 1.0, 1.0}}; }
    
    o->logScale = argparser_keyword_present(pp, "-logScale");
    
    o->reveal = argparser_keyword_present(pp, "-reveal");
    
    if (argparser_keyword_present(pp, "-shrink"))
      { o->shrink = argparser_get_next_double(pp, 0.0, 10.0); }
    else 
      { o->shrink = 1.0; }
    
    if (argparser_keyword_present(pp, "-exceptions"))
      { o->exceptions = argparser_get_next_non_keyword(pp); }
    else 
      { o->exceptions = NULL; }
    
    if (argparser_keyword_present(pp, "-stretched"))
      { o->stretched = argparser_get_next_non_keyword(pp); }
    else 
      { o->stretched = NULL; }
    
    int32_t k = 0; /* Layer count. */
    while (argparser_keyword_present(pp, "-layer"))
      { if (k >= LAYS) { argparser_error(pp, "too many layer specs"); }
        o->layer[k] = parse_layer_data(pp);
        k++;
      }
    if (k < LAYS) { argparser_error(pp, "too few layer specss"); }
    
    if (argparser_keyword_present(pp, "-bgColor"))
      { o->bgColor = frgb_parse_color(pp); }
    else 
      { /* The default background is layer 0: */
        o->bgColor = o->layer[0].color;
      }
    
    if (argparser_keyword_present(pp, "-outGamma"))
      { o->outGamma_expo = parse_triplet(pp, 0.1, 10.0); } 
    else 
      { o->outGamma_expo = (frgb_t){{1.0, 1.0, 1.0}}; }
    
    if (argparser_keyword_present(pp, "-outPrefix"))
      { o->outPrefix = argparser_get_next_non_keyword(pp); } 
    else 
      { o->outPrefix = "out/img"; }


    /* Get optional input file name: */
    if (argparser_next(pp) != NULL) 
      { o->f_in_name = argparser_get_next(pp); }
    else
      { o->f_in_name = "-"; }

    argparser_finish(pp);
        
    return o;
  }
  
layer_data_t parse_layer_data (argparser_t *pp)
  {
    layer_data_t dt;
    
    dt.name = argparser_get_next_non_keyword(pp);
    /* Parse layer's color (mandatory) */
    argparser_get_keyword_next(pp, "-color");
    dt.color = frgb_parse_color(pp);
      
    /* Parse new color for color stretching (optional): */
    if (argparser_keyword_present_next(pp, "-remap"))
      { dt.remap = frgb_parse_color(pp); }
    else
      { dt.remap = NoColor; }

    return dt;
  }
  
frgb_t parse_triplet(argparser_t *pp, double lo, double hi)
  { frgb_t p;
    for (int32_t i = 0;  i < 3; i++)
      { p.c[i] = (float)argparser_get_next_double(pp, lo, hi); }
    return p;
  }

void debug_color(char *label, int32_t col, int32_t row, frgb_t *fv, char *tail)
  { 
    if (main_debug) 
      { if (row >= 0) 
          { fprintf(stderr, "%s[%3d][%3d] = ", label, row, col); }
        else
          { fprintf(stderr, "%s = ", label); }
        print_color(stderr, "( ", fv, " )");
        fprintf(stderr, "%s", tail);
      }
  }

