#define PROG_NAME "pgmcolorize"
#define PROG_DESC "convert a grayscale image to pseudo-color"
#define PROG_VERS "1.0"

/* Copyright © 2006 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2024-12-21 12:00:53 by stolfi */

#define MAX_UNSIGNED_STYLE 0
#define MAX_SIGNED_STYLE 2
#define MAX_CYCLES 1000000
#define OUTPUT_MAXVAL PNM_FILE_MAX_MAXVAL

#define BT_BIAS (sample_conv_gamma_BT709_BIAS)

#define PROG_HELP \
  PROG_NAME "\\\n" \
  "    [ -zero {ZERO_REF} [ / {ZERO_SCALE} ] ] \\\n" \
  "    [ -inGamma {INGAMMA} ] \\\n" \
  "    [ -outGamma {OUTGAMMA} ] \\\n" \
  "    [ -cycles {NCYCLES} ] \\\n" \
  "    [ -style {STYLE} ] \\\n" \
  "    [ -showMap ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    < {INFILE} \\\n" \
  "    > {OUTFILE}"
  
#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads a grayscale image file and writes a colorized" \
  " version of it, where different gray values have been" \
  " replaced by different colors.  The colors mapping" \
  " is chosen so as to enhance small changes and help human" \
  " viewers to identify regions with similar values, independently" \
  " of their position and surrounds.\n" \
  "\n" \
  "FILE FORMATS\n" \
  "  The input image must be monochromatic, in PBM or PGM format.  The output" \
  " is always in PPM format with {MAXVAL = " stringify(OUTPUT_MAXVAL) "}\n" \
  "\n" \
  "INPUT SCALING\n" \
  "  Input image samples are first linearly scaled from [0..{MAXVAL_IN}]" \
  " to [0_1],where {MAXVAL_IN} is the nominal maximum pixel value declared" \
  " in the file's header.  Then, if the option \"-zero {Z}\" is specified," \
  " they are mapped again from the range [0_1] to the range [-1_+1], by" \
  " a piecewise linear map.  More precisely, the sub-range [0_{Z}] gets" \
  " mapped to [-1_0], and [{Z}_1] gets mapped to [0_+1], both linearly.\n" \
  "\n" \
  "COLOR MAPPING\n" \
  "  In either case, the resulting fractional" \
  " values are mapped to RGB color values, along a spiral path in the" \
  " RGB cube.  The color map is such that the luminance {Y} of the colors strictly" \
  " increases along the path, while cycling through several hues.\n" \
  "  If the \"-zero\" option is given, input values in [-1_0] are mapped to cold" \
  " hues, zero is mapped to gray, and values in [0_+1] are mapped to warm" \
  " hues.  Otherwise, both cold and warm hues are used for the range [0_1].\n" \
  "\n" \
  "GAMMA CORRECTION\n" \
  "  If the \"-inGamma\" option is specified, the scaled value {V} of each" \
  " pixel gets further mapped by the nonlinear decoding function" \
  " {sample_conv_gamma}, with decoding exponent {gammaDec=INGAMMA} and a" \
  " standard {bias} parameters.  Conversely, if \"-outGamma\" is" \
  " specified, each RGB coordinate {V} of each output pixel gets" \
  " mapped by {sample_conv_gamma} with enconding exponent {gammaEnc=1/OUTGAMMA}, before" \
  " quantization.  A typical value for {INGAMMA} and {OUGAMMA} (as" \
  " per ITU-R BT.709) is {1/0.450 = 2.222...}.\n" \
  "\n" \
  "  " sample_conv_gamma_INFO "\n" \
  "OPTIONS\n" \
  "\n" \
  "  -zero {ZERO_REF} [ / {ZERO_SCALE} ]\n" \
  "    Specifies that the input image value {Z = ZERO_REF/ZERO_SCALE}" \
  " should be mapped to zero.  Also specifies" \
  " that cold and warm hues should be used to denote negative and" \
  " positive values, respectively.  The default for {Z} is zero," \
  " meaning that the input values are interpreted as" \
  " non-negative.  If {ZERO_SCALE} is not specified, it defaults to 1." \
  " defaults to 1.\n" \
  "\n" \
  "  -style {STYLE}\n" \
  "    Specifies the style of color map to use. The argument {STYLE}" \
  " is a number between 0 and " stringify(MAX_UNSIGNED_STYLE) " for" \
  " unsigned input, and between 0 and " stringify(MAX_SIGNED_STYLE) " for" \
  " signed input.  The default is 0.\n" \
  "\n" \
  "  -cycles {NCYCLES}\n" \
  "    Specifies how many times the hues of the color map should cycle" \
  " through the spectrum. A negative value causes the hues to cycle" \
  " in the opposite sense  Note that this parameter must be 1 for" \
  " signed color scales (with positive {ZERO_REF}) and" \
  " {STYLE} 0.  The default is 1. \n" \
  "\n" \
  "  -inGamma {INGAMMA} \n" \
  "    Specifies the exponent for gamma-decoding the input" \
  " image file's samples.  The default is 1 (linear encoding).\n" \
  "\n" \
  "  -outGamma {OUTGAMMA} \n" \
  "    Specifies the exponent for gamma-decoding the output" \
  " image file's samples.  The default is 1 (linear encoding).\n" \
  "\n" \
  "  -showMap\n" \
  "    Requests a printout of the color map to {stderr}.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  ppmtopgm(1), pnmgamma(1), frgb_path.h.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created nov/2006 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP)."

#define stringify(x) strngf(x)
#define strngf(x) #x

#include <values.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <argparser.h>
#include <bool.h>
#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <frgb_path.h>
#include <sample_conv.h>

typedef struct ppm_pixel_t { uint16_t c[3]; } ppm_pixel_t;
  /* A RGB pixel value. */

typedef struct pnm_header_t 
  { pnm_format_t format;
    int32_t rows;
    int32_t cols;
    int32_t chns;
    uint16_t maxval;
    bool_t raw;          /* TRUE iff raw variant. */
    bool_t bits;         /* TRUE iff PBM format. */
  } pnm_header_t;
  /* Data for/from a PBM/PPM/PGM file header. */

typedef struct options_t 
  { char *filename_in;   /* Input filename, or "-" for stdin. */
    double inGamma_expo;      /* Nonlinearity exponent (gamma) of input file. */
    double outGamma_expo;     /* Nonlinearity exponent (gamma) of output file. */
    double zero;         /* Input value that corresponds to zero level. */
    int32_t style;           /* Style of color map to use. */
    int32_t cycles;          /* How many cycles through the spectrum. */
    bool_t showMap;      /* TRUE prints the colormap to {stderr}. */
  } options_t;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);

options_t *get_options(int32_t argc, char **argv);

void select_colors
  ( ppm_pixel_t *map, 
    int32_t maxval, 
    double zero, 
    int32_t style, 
    int32_t cycles, 
    double inGamma_expo, 
    double outGamma_expo,
    bool_t showMap
  );
  /* Stores in {map[0..nmap-1]} a suitable palette of RGB colors
    for input values {0..maxval}.
    
    If {zero <= 0}, assumes that the input pixel values in {0..maxval}
    represent fractions in the range {[0_1]}. Then, {map[0]} is black,
    {map[nmap-1]} is white, and intermediate entries have intermediate
    intensities (with varying hues).
    
    If {zero > 0}, assumes that the input pixels values from
    {0..maxval} represent fractions in the range {[-1_+1]}, with value
    {zero*maxval} corresponding to 0. Then {map[k]} is a cold color
    for {k < zero*maxval}, a warm color for {k > zero*maxval}, and
    gray for {k = zero*maxval}. In any case the luminosities increase
    with {k}.
    
    The {style} parameter chooses between several styles of color map,
    which differ in subtle details.
    
    In each range (positive or negative) the colors will go through a
    set of hues, {cycles} times. If {cycles} is zero, the colors will
    lie along a straight line in RGB space.
    
    The input values and the output RGB coordinates are adjusted
    according to the given {inGamma_expo} and {outGamma_expo}, respectively. If
    {zero > 0}, the input values are gamma-adjusted *after* mapping to
    the range [-1_+1]. */

void read_input_file_header(FILE *f, pnm_header_t *hd);
  /* Reads the header of a PBM/PGM/PPM file, returns the data in {hd}. */
  
void write_output_file_header(FILE *f, pnm_header_t *hd, pnm_header_t *hd_in);
  /* Creates the header {hd} of the output file, given that {hd_in} 
    of the input file, and writes it to file {f}. */

void process_pixels(FILE *f, pnm_header_t *hd_in, options_t *o, ppm_pixel_t *map, pnm_header_t *hd_ot);
  /* Reads body of input image from file {f}, conerts to grayscale,
    applies the pseudocolor map {map}, and writes the color pixels to
    stdout.  Adjust for input and output gamma. */

double gamma_decode(double y, double gammaDec);
  /* Returns {sample_conv(y,gammaDec,BT_BIAS)}, which is 
    roughly {sign(y)*abs(y)^gammaDec}. */
  
double gamma_encode(double y, double gammaDec);
  /* Returns {sample_conv(y,1/gammaDec,BT_BIAS)}, which is 
    roughly {sign(y)*abs(y)^(1/gammaDec)}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    pnm_header_t hd_in;           /* Specs from input file header. */
    pnm_header_t hd_ot;           /* Specs for output file header. */

    /* Open input file and read header: */
    FILE *rd = ((o->filename_in == NULL) ? stdin : open_read(o->filename_in, TRUE));
    read_input_file_header(rd, &hd_in);

    /* Open output file, choose output {maxval} and {format} and write header: */
    /* Choose maxval and format for output image: */
    fprintf(stderr, "writing result to (stdout)...\n");
    write_output_file_header(stdout, &hd_ot, &hd_in);

    /* Allocate and set up the color table {map[0..ncolors-1]}: */
    int32_t ncolors = hd_in.maxval + 1;
    ppm_pixel_t *map = (ppm_pixel_t *)notnull(malloc(ncolors*sizeof(ppm_pixel_t)), "no mem");
    select_colors
      ( map, hd_in.maxval, 
        o->zero, o->style, o->cycles, 
        o->inGamma_expo, o->outGamma_expo, 
        o->showMap
      );
    
    /* Apply color table to all pixels: */
    process_pixels(rd, &hd_in, o, map, &hd_ot);

    fflush(stdout);
    if (rd != stdin) { fclose(rd); }
    fprintf(stderr, "done.\n");
    return 0;
  }
  
void process_pixels
  ( FILE *f, 
    pnm_header_t *hd_in, 
    options_t *o,
    ppm_pixel_t *map, 
    pnm_header_t *hd_ot
  )
  {
    /* File format consistency checks: */
    int32_t fmt_in = hd_in->format;
    if 
      ( (fmt_in != PBM_FORMAT) && (fmt_in != RPBM_FORMAT) && 
        (fmt_in != PGM_FORMAT) && (fmt_in != RPGM_FORMAT)
      )
      { pnm_error("input image must be monochromatic (PBM or PGM)"); exit(1); }
    
    int32_t fmt_ot = hd_ot->format;
    if ((fmt_ot != RPPM_FORMAT) && (fmt_ot != PPM_FORMAT))
      { pnm_error("output image must be color (PPM)"); exit(1); }

    int32_t rows = hd_in->rows;
    int32_t cols = hd_in->cols;
    assert(hd_in->chns == 1);
    assert(hd_ot->chns == 3);
    
    /* Read pixels: */
    uint16_t *smp_in = uint16_image_alloc_pixel_row(cols, hd_in->chns);
    uint16_t *smp_ot = uint16_image_alloc_pixel_row(cols, hd_ot->chns);
    int32_t row, col;
    for (row = 0; row < rows; ++row)
      { pnm_read_pixels(f, smp_in, cols, hd_in->chns, hd_in->maxval, hd_in->raw, hd_in->bits);
        uint16_t *p_in = smp_in;
        uint16_t *p_ot = smp_ot;
        for (col = 0; col < cols; ++col)
          { /* Extract pixel {iv} from input PGM/PBM image: */
            int32_t iv = *p_in; ++p_in;
            /* if (row == 0) { fprintf(stderr, "row[%4d] = %5d\n", col, iv); } */
            /* Convert to pseudocolor: */
            ppm_pixel_t px = map[iv];
            /* Save in output PPM image: */
            (*p_ot) = px.c[0]; p_ot++;
            (*p_ot) = px.c[1]; p_ot++;
            (*p_ot) = px.c[2]; p_ot++;
          }
        pnm_write_pixels(stdout, smp_ot, cols, hd_ot->chns, hd_ot->maxval, hd_ot->raw, hd_ot->bits);
      }
  }

void read_input_file_header(FILE *f, pnm_header_t *hd)
  { pnm_read_header
      ( f, &(hd->cols), &(hd->rows), &(hd->chns), 
        &(hd->maxval), &(hd->raw), &(hd->bits), 
        &(hd->format)
      );
    fprintf(stderr, "input format = %i%i maxval = %d cols = %d rows = %d\n", 
      (hd->format >> 8)&255, hd->format&255, hd->maxval, hd->rows, hd->cols);
  }

void write_output_file_header(FILE *f, pnm_header_t *hd, pnm_header_t *hd_in)
  { hd->cols = hd_in->cols;
    hd->rows = hd_in->rows;
    hd->chns = 3;
    hd->maxval = OUTPUT_MAXVAL;
    hd->raw = TRUE;
    hd->bits = FALSE;
    hd->format = RPPM_FORMAT;
    pnm_write_header(f, hd->cols, hd->rows, hd->maxval, hd->format);
  }

void select_colors
  ( ppm_pixel_t *map, 
    int32_t maxval, 
    double zero,
    int32_t style, 
    int32_t cycles, 
    double inGamma_expo,
    double outGamma_expo,
    bool_t showMap
  )
  {
    double scale = (double)(maxval);
    int32_t print_step = (maxval <= 255 ? 1 : (maxval + 1)/256); 
    int32_t vi;
    for (vi = 0; vi <= maxval; vi++)
      { /* Compute float value {fi} of pixel value {vi}, in linear [0_1] range: */
        double fi = vi/scale;
        
        /* Convert {fi} to {zi} in desired range ([0_1] or [-1_+1]): */
        double zi = fi - zero;
        if (zero <= 0) 
          { /* Unsigned input -- rescale to keep in [0_1]: */
            zi = zi/(1 - zero);
            assert(zi >= 0);
            assert(zi <= 1);
          }
        else
          { /* Signed input -- remap to keep in [-1_+1]: */
            zi = zi/(zi <= 0 ? zero : 1 - zero);
            assert(zi >= -1);
            assert(zi <= +1);
          }

        /* Apply input gamma correction, result is {gi}: */
        double gi = gamma_decode(zi, inGamma_expo);

        /* Compute output color {go} as linear float RGB: */
        frgb_t go;
        if (zero <= 0)
          { /* Color maps for unsigned values: */
            switch(style)
              { case 0: go = frgb_path_map_unsigned_0(gi, cycles); break;
                default: assert(FALSE); 
              }
          }
        else
          { /* Color maps for signed values: */
            switch(style)
              { case 0: 
                  if (cycles != 1) 
                    { fprintf(stderr, "{cycles} = %d ignored for signed path type 0", cycles); }
                  go = frgb_path_map_signed_0(gi, cycles); break;
                case 1:
                  go = frgb_path_map_signed_1(gi, cycles); break;
                case 2: 
                  go = frgb_path_map_signed_2(gi, cycles); break;
                default: 
                  assert(FALSE); 
              }
          }
          
        /* Apply output gamma to obtain nonlinear float RGB color {zo}: */
        frgb_t zo;
        zo.c[0] = (float)gamma_encode(go.c[0], outGamma_expo);
        zo.c[1] = (float)gamma_encode(go.c[1], outGamma_expo);
        zo.c[2] = (float)gamma_encode(go.c[2], outGamma_expo);
        
        /* Quantize coords as {vo.c[0..2]}: */
        ppm_pixel_t vo;
        bool_t isMask = FALSE; /* Quantize as a smooth-valued image. */
        vo.c[0] = pnm_quantize(zo.c[0], OUTPUT_MAXVAL, isMask, PNM_NO_BADVAL);
        vo.c[1] = pnm_quantize(zo.c[1], OUTPUT_MAXVAL, isMask, PNM_NO_BADVAL);
        vo.c[2] = pnm_quantize(zo.c[2], OUTPUT_MAXVAL, isMask, PNM_NO_BADVAL);
        
        /* Store in map: */
        map[vi] = vo;
        
        /* Print color map: */
        bool_t interesting = ((vi == maxval) || (vi % print_step == 0));
        if (showMap && interesting)
          { fprintf(stderr, "%5d = %6.4f = %+7.4f = %+7.4f", vi, fi, zi, gi);
            fprintf(stderr, " --> ( %5.3f %5.3f %5.3f )", go.c[0], go.c[1], go.c[2]);
            fprintf(stderr, " [ %5.3f ]", frgb_get_Y(&go));
            fprintf(stderr, " --> ( %5.3f %5.3f %5.3f )", zo.c[0], zo.c[1], zo.c[2]);
            fprintf(stderr, " --> ( %03d %03d %03d )", vo.c[0], vo.c[1], vo.c[2]);
            fprintf(stderr, "\n");
          }
      }
  }

double gamma_decode(double y, double gammaDec)
  {
    return sample_conv_gamma((float)y, gammaDec, BT_BIAS);
  }

double gamma_encode(double y, double gammaDec)
  {
    return sample_conv_gamma((float)y, 1/gammaDec, BT_BIAS);
  }

options_t *get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    if (argparser_keyword_present(pp, "-zero"))
      { o->zero = argparser_get_next_double(pp, -1.0e100, 1.0e100);
        if (argparser_keyword_present_next(pp, "/"))
          { double den = argparser_get_next_double(pp, 1.0e-100, 1.0e100);
            o->zero /= den; 
          }
      }
    else 
      { o->zero = 0.0; }
    
    if (argparser_keyword_present(pp, "-style"))
      { int32_t max_style;
        if (o->zero <= 0)
          { max_style = frgb_path_map_unsigned_max_style(); }
        else
          { max_style = frgb_path_map_signed_max_style(); }
        o->style = (int32_t)argparser_get_next_int(pp, 0, max_style);
      }
    else 
      { o->style = 0; }
    
    if (argparser_keyword_present(pp, "-cycles"))
      { o->cycles = (int32_t)argparser_get_next_int(pp, -(MAX_CYCLES), +(MAX_CYCLES)); }
    else 
      { o->cycles = 1; }
      
    /* Some path types do not support multiple cycles: */
    if (o->cycles != 1)
      { if (o->zero <= 0)
          { if (o->style == 1)
              { argparser_error(pp, "{cycles} must be 1 for unsigned input, style 1"); }
          }
        else
          { if (o->style == 0)
              { argparser_error(pp, "{cycles} must be 1 for signed input, style 0"); }
          }
      }
         
    
    if (argparser_keyword_present(pp, "-inGamma"))
      { o->inGamma_expo = argparser_get_next_double(pp, 1.0e-5, 1.0e5); }
    else 
      { o->inGamma_expo = 1.0; }
    
    if (argparser_keyword_present(pp, "-outGamma"))
      { o->outGamma_expo = argparser_get_next_double(pp, 1.0e-5, 1.0e5); }
    else 
      { o->outGamma_expo = 1.0; }
      
    o->showMap = argparser_keyword_present(pp, "-showMap");
    
    argparser_finish(pp);

    return o;
  }
  
/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 2004 by the State University of Campinas (UNICAMP).
**
** Created in 2004 by Jorge Stolfi, IC-UNICAMP.
** Based on color map function from jclimage.c by Jorge Stolfi.
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
