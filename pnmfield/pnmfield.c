#define PROG_NAME "pnmfield"
#define PROG_DESC "generate a PGM/PPM image with a wave-like color field"
#define PROG_VERS "1.0"

#define pnmfield_C_COPYRIGHT "Copyright © 2003 by the State University of Campinas (UNICAMP)"

/* Last edited on 2024-12-20 18:05:28 by stolfi */

/* TO DO: !!! Unify the documentation with that of {pnmadjust.c} !!! */

#define PROG_HELP \
  PROG_NAME "\\\n" \
  "  [ -inGamma Nr Ng Nb ].. \\\n" \
  "  -field \\\n" \
  "    { uniform         RGB0 \\\n" \
  "      hWave       H0  RGB0   H1   RGB1 | \\\n" \
  "      vWave       V0  RGB0   V1   RGB1 | \\\n" \
  "      wave      H0 V0 RGB0  H1 V1 RGB1 | \\\n" \
  "      wavePair  H0 V0 RGB0  H1 V1 RGB1  H2 V2 RGB2 \\\n" \
  "      ramp      H0 V0 RGB0  H1 V1 RGB1  H2 V2 RGB2 \\\n" \
  "    } \\\n" \
  "  [ -logarithmic ].. \\\n" \
  "  [ -gray ].. \\\n" \
  "  [ -maxval PIXVAL ].. \\\n" \
  "  [ -outGamma Nr Ng Nb ].. \\\n" \
  "  WIDTH HEIGHT\n" \
  "where the Hi, Vi are pixel indices, each Ni is a positive coefficient,\n" \
  "and each RGBi is three intensity values Ir Ig Ib, optionally followed by\n" \
  "a slash \"/\" and a common denominator Iw."

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Writes a PPM image of {WIDTH} columns and {HEIGHT} rows" \
  " that contains a smoothly-varying color field.\n" \
  "\n" \
  PROG_INFO_OPTIONS "\n" \
  "\n" \
  PROG_INFO_FIELD_SPECS "\n" \
  "\n" \
  PROG_INFO_EXAMPLE "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmgamma(1), ppmxramp(1), ppmmake(1), pgmnorm(1), pnm(5)\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created apr/2003 by Jorge Stolfi, IC-UNICAMP,\n" \
  " based on ColorCorrect.m3 by Jorge Stolfi, jun/1995.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2008-07-07 by J. Stolfi, IC-UNICAMP.  Folded th emanpage into the program.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmfield_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_OPTIONS \
  "OPTIONS\n" \
  "  In what follows, the {TRIPLET} parameters are three" \
  " real numbers, one for each color channel (red, green," \
  " blue).  The {COLOR} parameters are " frgb_parse_color_INFO "\n" \
  "\n" \
  "  The {COLOR} components refer to sample values in the" \
  " input image, only scaled to the [0 _ 1] range; they are thus assumed to" \
  " be affected by light-field factors, dark-level noise, and gamma" \
  " encoding.\n" \
  "\n" \
  "  -inGamma {TRIPLET}\n" \
  "    Specifies the exponent {expo} to be assumed when" \
  " interpreting user-given color values. Namely, each" \
  " component of the color is converted to a number {y}" \
  " between 0 and 1, and then mapped through" \
  " {sample_conv_gamma(y,expo,bias)} where {expo}" \
  " is the corresponding component of the {TRIPLET} and {bias} is" \
  " {sample_conv_BT709_BIAS}.  (For typical monitors, {expo}" \
  " usually ranges between 1.8 and 2.2.)  The default is 1.0.\n" \
  "\n" \
  "  -field {FIELDSPECS}\n" \
  "    Specifies the kind of color field and its parameters;" \
  " see the FIELD SPECS section below.\n" \
  "\n" \
  "  -logarithmic\n" \
  "    Specifies that the waves should be senoidal in the" \
  " logarithmic scale.\n" \
  "\n" \
  "  -outGamma {TRIPLET}\n" \
  "    Specifies the nonlinear exponent to be used when encoding" \
  " the corrected and adjusted intensities before writing them" \
  " to the output file. See the \"-inGamma\" option for details.\n" \
  "\n" \
  "  -gray\n" \
  "    Specifies PGM (grayscale) output instead of PPM. (Note that" \
  " all color and gamma parameters must be given with three RGB" \
  " components, even in this case.)\n" \
  "\n" \
  "  -maxval {NUMBER}\n" \
  "    The maximum pixel value for each component of the output" \
  " image.  The default is 255 for color images, and" \
  " {2^16-1 = 65535} for   grayscale."

#define PROG_INFO_FIELD_SPECS \
  "FIELD SPECS\n" \
  "  The arguments following the \"-fields\" keyword" \
  " may have the following format:\n" \
  "\n" \
  "    uniform {COLOR}\n" \
  "      Specifies a uniform image of the given color.\n" \
  "\n" \
  "    ramp {H0} {V0} {COLOR0}  {H1} V1} {COLOR1} {H2} {V2} {COLOR2}\n" \
  "      Specifies a field that varies linearly with position," \
  " and has value {COLORi} at column {Hi} and row {Vi}, for  {i = 0,1,2}.\n" \
  "\n" \
  "    hWave {H0} {COLOR}0 {H1} {COLOR1}\n" \
  "      Specifies a simple senoidal wave travelling" \
  " horizontally.  Argument {H0} is a pixel column where" \
  " the wave is maximum, and {H1} the nearest coordinate" \
  " where the wave is minimum.  The colors {COLOR0} and" \
  " {COLOR1} are the corresponding values.\n" \
  "\n" \
  "    vWave {V0} {COLOR0} {V1} {COLOR1}\n" \
  "      Analogous to \"hWave\" except that the wave travels" \
  " vertically. {V0} and {V1} are the row indices of consecutive" \
  " wave maximum and minimum.\n" \
  "\n" \
  "    wave {H0 V0 COLOR0 {H1 V1 COLOR1}\n" \
  "      Specifies a single senoidal wave in a general" \
  " direction.  Arguments {H0} {V0} are the coordinates of a" \
  " pixel where the wave is maximum, {H1} {V1} is the nearest" \
  " pixel where the wave is minimum, and {COLOR0}, {COLOR1} are" \
  " the corresponding color values.\n" \
  "\n" \
  "    wavePair {H0} {V0} {COLOR0}  {H1} {V1} {COLOR1} {H2} {V2} {COLOR2}\n" \
  "      Specifies the combination of two senoidal waves in general" \
  " directions. Arguments {H0} {V0} are the coordinates of a pixel" \
  " where the both waves are maximum, {H1} {V1} is the nearest pixel" \
  " where wave 1 is at its minimum, and {H2} {V2} is the nearest pixel" \
  " where wave 2 is at its minimum.  The image will have value {COLOR0}" \
  " applies where both waves are at a maximum, {COLOR1} where wave 1 is" \
  " maximum and wave 2 is minimum, and {COLOR2} where wave 2 is maximum" \
  " and wave 1 is  minimum.  (Note that these last two parameters are" \
  " the colors of pixels {H1} {V1} and {H2} {V2} only if the waves are" \
  " perpendicular to each other.)\n" \
  "\n" \
  "  In all these specifications, the pixel coordinates  {Hi}, {Vi}" \
  " may lie outside the image."

#define PROG_INFO_EXAMPLE \
  "EXAMPLE\n" \
  "  The following example generates a 640 by 480 image that" \
  " varies according to two orthogonal waves, with periods" \
  " (300,100) and (-50,150):\n" \
  "\n" \
  "    pnmfield \\\n" \
  "      -inGamma      1.0 1.0 1.0       \\\n" \
  "      -field wavePair  \\\n" \
  "         0120 0320  240 235 230 / 255 \\\n" \
  "         0420 0420  201 198 180 / 255 \\\n" \
  "         0070 0470  185 180 178 / 255 \\\n" \
  "         640 480"

#include <values.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include <argparser.h>

#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>

#include <frgb.h>
#include <frgb_ops.h>
#include <colorfield.h>
#include <sample_conv.h>

/* Maximum rows, columns, elements (to avoid absurd allocs): */
#define MAX_IMG_PIXELS (16*256*256*256)   
#define MAX_IMG_SIZE MAX_IMG_PIXELS

#define MAXCHANNELS 3

#define BT_BIAS (sample_conv_BT709_BIAS) 
/* A {bias} parameter that approximates BT.709 when
  used with expo {0.45} or {1/0.45}. */


/* COMMAND LINE ARGUMENTS */

typedef struct options_t 
  { frgb_t inGamma_expo;      /* Gamma of given color values, per channel. */
    int cols, rows;      /* Dimensions of output image. */
    cfld_args_t *fld;    /* Raw user specs for the color field. */
    bool_t logarithmic;  /* TRUE to compute the field in log scale. */
    bool_t gray;         /* Write a grayscale image. */
    frgb_t outGamma_expo;     /* Gamma of output image, per channel. */
    uint16_t maxval; /* Max pixel value on output. */
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments. */
  
void write_pixels(FILE *f, options_t *o);
  /* Computes pixels, writes them to file {f}. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  { options_t *o = parse_options(argc, argv);
    write_pixels(stdout, o);
    fflush(stdout);
    return 0;
  }
   
void write_pixels(FILE *f, options_t *o)
  {
    uint16_t maxval = o->maxval;
    int rows = o->rows;
    int cols = o->cols;
    int chns = (o->gray ? 1 : 3);
    pnm_format_t format;
    bool_t raw, bits;
    pnm_choose_output_format(o->maxval, chns, FALSE, &format, &raw, &bits);
    pnm_write_header(f, cols, rows, maxval, format);

    int needs_inGamma = ! frgb_fequal(o->inGamma_expo.c, frgb_Ones.c, chns);
    frgb_t *inGamma_expo = (needs_inGamma ? &(o->inGamma_expo) : NULL);

    int needs_outGamma = ! frgb_fequal(o->outGamma_expo.c, frgb_Ones.c, chns);
    frgb_t *outGamma_expo = (needs_outGamma ? &(o->outGamma_expo) : NULL);

    auto frgb_t adjust_arg(frgb_t *v, int col, int row);
    
    frgb_t adjust_arg(frgb_t *v, int col, int row)
      { return frgb_correct_arg(v, inGamma_expo, o->gray); }
    
    cfld_params_t *fld = cfld_compute_params(o->fld, adjust_arg, o->logarithmic);

    /* Input and output scaling parameters: */
    double zero = 0, unit = (double)maxval;
    double scale = unit - zero;
    
    uint16_t *smp = uint16_image_alloc_pixel_row(cols, chns);
    uint16_t *sP;
    
    int iv[MAXCHANNELS];
    int y, x, c;
    for (y = 0; y < rows; ++y)
      { for(x = 0, sP = smp; x < cols;  ++x)
          { frgb_t fv;
            frgb_DEBUG = (x == 0) & (y == 0);
            cfld_eval(fld, x, y, &fv, chns);
            if (chns == 3)
              { frgb_clip_rgb(&fv); }
            else
              { fv.c[0] = frgb_clip_gray(fv.c[0]); }
            frgb_debug("fv", x, y, &fv, chns, "\n");
            if (outGamma_expo != NULL)
              { for (c = 0; c < chns; c++)
                  { fv.c[c] = sample_conv_gamma(fv.c[c], 1/outGamma_expo->c[c], BT_BIAS); }
              }
            frgb_debug("op", x, y, &fv, chns, "\n");
            for (c = 0; c < chns; c++) 
              { iv[c] = frgb_quantize(fv.c[c], zero, scale, maxval); 
                (*sP) = (uint16_t)iv[c]; 
                sP++;
              }
          
            frgb_debug_int_pixel("oq", x, y, iv, chns, "\n");
            if (frgb_DEBUG) { fprintf(stderr, "\n\n"); }
          }
        pnm_write_pixels(f, smp, cols, chns, maxval, raw, bits);
      }
  }

options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
     
    options_t *o = (options_t *)malloc(sizeof(options_t));
    
    /* Parse keyword arguments: */
    if (argparser_keyword_present(pp, "-inGamma"))
      { o->inGamma_expo = frgb_parse(pp, 0.1, 10.0); }
    else 
      { o->inGamma_expo = (frgb_t){{1.0, 1.0, 1.0}}; }
   
   if (argparser_keyword_present(pp, "-field"))
      { o->fld = cfld_parse(pp); }
    else 
      { o->fld = NULL; }
      
    o->gray = argparser_keyword_present(pp, "-gray");
    
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = argparser_get_next_int(pp, 0, PNM_FILE_MAX_MAXVAL); }
    else 
      { o->maxval = (o->gray ? 65535 : 255 ); }
    
    
    o->logarithmic = argparser_keyword_present(pp, "-logarithmic");
    
    if (argparser_keyword_present(pp, "-outGamma"))
      { o->outGamma_expo = frgb_parse(pp, 0.1, 10.0); } 
    else 
      { o->outGamma_expo = (frgb_t){{1.0, 1.0, 1.0}}; } 

    /* Parse X and Y sizes: */
    argparser_skip_parsed(pp);
    o->cols = argparser_get_next_int(pp, 0, MAX_IMG_SIZE);
    o->rows = argparser_get_next_int(pp, 0, MAX_IMG_SIZE);
    if (MAX_IMG_PIXELS/o->cols < o->rows) 
      { pnm_error("too many pixels"); }
    
    /* Check for extraneous arguments: */
    argparser_finish(pp);
        
    return o;
  }
