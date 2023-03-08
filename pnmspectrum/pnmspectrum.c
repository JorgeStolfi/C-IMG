#define PROG_NAME "pnmspectrum"
#define PROG_DESC "computes the power spectrum of a PBM/PGM/PPM image"
#define PROG_VERS "1.0"

/* Last edited on 2023-01-14 12:09:24 by stolfi */

#define pnmspectrum_C_COPYRIGHT \
  "Copyright © 2008 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    [ -outputImage { {PNM_NAME} | NONE } ] \\\n" \
  "      [ -center ] \\\n" \
  "      [ -symmetric ] \\\n" \
  "      [ -scale [ linear | log {MIN_PWR} ] ] \\\n" \
  "      [ -maxPower { {MAX_PWR} | AUTO } ] \\\n" \
  "      [ -maxPixel {MAX_PIX} ] \\\n" \
  "    [ -outputTable { {TXT_NAME} | NONE } ] \\\n" \
  "      [ -exact | -freqRanges {NFREQ} ] \\\n" \
  "    [ -zeroMean ] \\\n" \
  "    [ -averageImages ] \\\n" \
  "    [ -vignette ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    [ {PNMFILE_IN} .. ]"

#define PROG_INFO_DESC \
  "  The program computes the average Fourier (actually, Hartley)" \
  " power spectrum of one or more PNM images {PNMFILE_IN}, all with the same width and height.  The" \
  " result can be written out as a two-dimensional PNM image, isologous to" \
  " the Hartley transform, or as a text file containg a one-dimensional" \
  " table, that gives the power as a function of absolute spatial frequency.\n" \
  "\n" \
  "THE HARTLEY TRANSFORM\n" \
  float_image_hartley_INFO "\n" \
  "\n" \
  "THE HARTLEY POWER SPECTRUM\n" \
  "  The Hartley power spectrum is simply the Hartley transform with each" \
  " sample squared.  Thus, the sample in column {FX} and row {FY} is the" \
  " power contained in the Hartley component with those" \
  " frequencies.  The \"-center\" and \"-symmetric\" options" \
  " may make the spectrum image easier to understand."

#define PROG_INFO_DESC_IMG_OUT \
  "   The output spectrum image has the same dimensions {NX,NY} and the same number" \
  " of channels {NC} as the original image.  For each channel, the pixel in" \
  " column {FX} (from left) and row {FY} (from bottom) is the" \
  " square of the Hartley transform component with horizontal" \
  " frequency {FX} and vertical frequency {FY}, averaged over all images.  The output" \
  " image file is in PGM format if {NC==1}, and in PPM format if {NC==3}.\n" \
  "\n" \
  "  Each sample {V_OUT} of the spectrum is rescaled and quantized" \
  " by {sample_conv_quantize(V_OUT,MV_OUT,FALSE,V_LO,V_HI,...)}."

#define PROG_INFO_DESC_TXT_OUT \
  "   The TXT output contains a few header" \
  " lines, beginning with '#', and then one or more data" \
  " lines.  Each data line contains {} fields\n" \
  "\n" \
  "   \"{FMIN} {FMAX}  {F}  {NTERMS} {POWER}\"\n" \
  "\n" \
  " where\n" \
  "\n" \
  "   {FMIN}, {FMAX} is a range of spatial frequencies (in waves per pixel),\n" \
  "\n" \
  "   {F} is some spatial frequency in that range,\n" \
  "\n" \
  "   {NTERMS} is the number of Hartley terms with" \
  " pixel frequencies in that range,\n" \
  "\n" \
  "   {POWER} is the sum of the squares of the coefficients of those" \
  " terms, over all selected channels.\n" \
  "\n" \
  "  The ranges {FMIN} {FMAX} are consecutive and increasing, i.e. the" \
  " {FMAX} of one line is equal to the {FMIN} of the next line.  The" \
  " first line has {FMIN = 0}.\n" \
  "\n" \
  "  The text output may have either of two formats, /exact/ or" \
  " /binned/, depending on whether the command line option \"-freqRanges\" is" \
  " given or not.  The exact variant is most informative, while the binned" \
  " variant is more suitable for plotting and for comparing the" \
  " spectra of different images (especially when they have" \
  " different dimensions).\n" \
  "\n"

#define PROG_INFO_DESC_TXT_OUT_EXACT \
  "  If the option \"-exact\" is specified, the output text file has one line" \
  " for each distinct pixel frequency {F = sqrt((FX/NX)^2+(FY/NY)^2)}" \
  " appearing in the Hartley transform; where any two frequencies that yield the" \
  " same value when converted to {float} are treated as equal.  These `proper' lines" \
  " have {FMIN == FMAX == F}.  The first line," \
  " for the constant term, has {FMIN == FMAX == F == 0}. Note that the values" \
  " of {F} are irregularly spaced.  The field {NTERMS} in each line (which is" \
  " always an integer, usually 4) is the" \
  " number of Hartley terms that have that frequency. There are also `filler'" \
  " lines with {NTERMS = NPOWER = 0}, inserted between and around the proper lines" \
  " so that the ranges {[FMIN _ FMAX]} are contiguous and the last line has\n" \
  " {FMAX == F == sqrt(0.5)}.  The {F} values in the filler lines are" \
  " approximately halfway between {FMIN} and {FMAX}."

#define PROG_INFO_DESC_TXT_OUT_BINNED \
  "  If the option \"-freqRanges {NFREQS}\" is specified, there" \
  " are exactly {NFREQS} lines, and the range breaks {FMIN,FMAX} are" \
  " chosen so that each line has approximately the same" \
  " value of {NTERMS}.  Furthermore, the count and power of each term is" \
  " assumed to be spread over the range of frequencies between" \
  " {FLO} and {FHI} where\n" \
  "   {FLO = sqrt(((FX - 1)/NX)^2+((FY - 1)/NY)^2)}\n" \
  "   {FHI = sqrt(((FY + 1)/NX)^2+((FY + 1)/NY)^2)}\n" \
  " but mostly near {F}.  In particular, the last line has\n" \
  "   {FMAX = sqrt((floor(NX/2 + 1)/NX)^2+(floor(NY/2 + 1)/NY)^2)}\n" \
  " which is slightly larger than {sqrt(0.5)}." \
  "\n" \
  "  It follows that {NTERMS} is usually fractional, and may be" \
  " nonzero even if no Hartley terms have pixel frequencies in" \
  " the range {FMIN _ FMAX}.  The value {F} is some arbitrary" \
  " value in that range.  When analyzing the spectrum of an" \
  " image, one may want to plot" \
  " the ratios {PWR[i]/NTERMS}, rather than {PWR[i]} alone."

#define PROG_INFO_DESC_IMG_IN \
  "  If the argument {PNMFILE_IN} is omitted or is \"-\", the" \
  " program reads the PBM/PGM/PPM input image from {stdin}.  Each" \
  " sample {V_IN} is converted to a floating-point value in {[0_1]}" \
  " by {sample_conv_floatize(V_IN,MV_IN,FALSE,0.0,1.0,...)}."

#define PROG_INFO_OPTS_1 \
  "  -outputTable {TXT_NAME}\n" \
  "  -outputTable NONE\n" \
  "    The first form of this option specifies that the radial" \
  " spectrum should be written out as a text file named \"{TBL_NAME}\".  The second" \
  " form supresses that output.  At least one of" \
  " the options \"-outputImage\" and \"-outputTable\" must" \
  " be specified.  The default is \"-outputTable NONE\".\n" \
  "\n" \
  "  -outputImage {PNM_NAME}\n" \
  "  -outputImage NONE\n" \
  "    The first form of this option specifies that the spectrum" \
  " should be written out as a PNM image file.  The output image" \
  " file name will be \"{PNM_NAME}\", and will be a PGM or PPM image" \
  " depending on whether the input images are monochromatic or color," \
  " respectively.  If {PNM_NAME} is \"-\", the spectrum is written to standard output. The second" \
  " form supresses the spectrum image output.  At least one of the" \
  " options \"-outputImage\" and \"-outputTable\" must" \
  " be specified.  The default is \"-outputImage NONE\"."

#define PROG_INFO_OPTS_2 \
  "\n" \
  "  -center\n" \
  "    If this option is given, the spectrum image is shifted" \
  " cyclically along both axes, so that the the constant" \
  " component, with frequency {0,0}, is stored in the center" \
  " pixel in column {floor(NX/2)} (from left) and" \
  " row {floor(NY/2)} (from bottom).\n" \
  "\n" \
  "  -symmetric\n" \
  "    If this option is given, the entries of frequency {FX,FX}" \
  "    and {-FX,-FY} (modulo {NX,NY}) with their arithmetic mean, thus" \
  "    making the spectrum image symmetric and erasing the phase" \
  "    information that is present in the raw spectrum.\n" \
  "\n" \
  "  -maxPower {MAX_PWR}\n" \
  "  -maxPower AUTO\n" \
  "    This option specifies {MAX_PWR} as the maximum value in the power" \
  " spectrum that will be converted to the maximum sample value in" \
  " the output image.  Any samples that" \
  " exceed {MAX_PWR} will be clipped to {MAX_PWR}.  If the \"AUTO\" variant" \
  " is used, or if {MAX_PWR} is zero, the program will set {MAX_PWR} to" \
  " the maximum sample value that actually occurs in the power spectrum" \
  " of all channels.  This option" \
  " does not affect the table output.  If this option" \
  " is not given, the program assumes \"-maxPower AUTO\".\n" \
  "\n" \
  "  -scale linear\n" \
  "  -scale log {MIN_PWR}\n" \
  "    This option controls the conversion of each value {P} in" \
  " the power spectrum to the correponding sample value in the" \
  " output image.  The \"linear\" variant requests linear conversion," \
  " meaning that the sample value will be proportional to {P}.  The" \
  " second variant requests logarithmic" \
  " conversion, meaning that the output image sample" \
  " will be proportional to  {log(P/MIN_PWR)}.  In this case," \
  " power values less than {MIN_PWR} will be mapped to 0, and" \
  " the maximum sample value for output image conversion will be" \
  " implicitly changed from {MAX_PWR} to" \
  " {log(MAX_PWR/MIN_PWR)}.  This option" \
  " does not affect the text output.  If this option" \
  " is not given, the program assumes \"-scale linear\".\n" \
  "\n" \
  "  -zeroMean\n" \
  "    If this option is given, the computed power spectrum is" \
  " adjusted by setting the element with zero frequency (the" \
  " image's mean value) to zero.  This change affects the scaling" \
  " of power values in the output (image or text) file.  If not" \
  " specified, the mean value is not supressed.\n" \
  "\n" \
  "  -avgerageImages\n" \
  "    If this option is given, the images are averaged before" \
  " computing the spectrum.  Otherwise the power spectrum is" \
  " computed for each image, and the spectra are averaged.\n" \
  "\n" \
  "  -vignette\n" \
  "    If this option is given, images are multiplied by" \
  " a two-dimensional Hann windowing function, before computing the" \
  " spectra.  This has the effect of blurring the spectrum" \
  " slightly but gets rid of the spikes due to discontinuities" \
  " at the edges.\n" \
  "\n" \
  "  -maxPixel {MAX_PIX}\n" \
  "    Specifies {MAX_PIX} as the maximum sample value for the" \
  " output image.  It must be an integer between 255 and 65535," \
  " inclusive. If not specified, it is set to 65535 (meaning 16-bit unsigned samples)."

#define PROG_INFO_OPTS_3 \
  "  -exact\n" \
  "    Specifies that the output table should have a separate" \
  " entry for each distinct frequency (in waves/pixel) occurring" \
  " in the Hartley spectrum.\n" \
  "\n" \
  "  -freqRanges {N_FREQS}\n" \
  "    Specifies that the output table should have {N_FREQS} bins," \
  " spaced so that each bin contains approximately the same number" \
  " of Hartley spectrum terms.\n" \
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
  "IMAGE OUTPUT\n" \
  PROG_INFO_DESC_IMG_OUT "\n" \
  "\n " \
  "TEXT OUTPUT\n" \
  PROG_INFO_DESC_TXT_OUT "\n" \
  "\n" \
  "Exact Text Output\n" \
  PROG_INFO_DESC_TXT_OUT_EXACT "\n" \
  "\n" \
  "Binned Text Output\n" \
  PROG_INFO_DESC_TXT_OUT_BINNED "\n" \
  "\n " \
  "IMAGE INPUT\n" \
  PROG_INFO_DESC_IMG_IN "\n " \
  "\n " \
  "OPTIONS\n" \
  PROG_INFO_OPTS_1 "\n" \
  "\n" \
  PROG_INFO_OPTS_2 "\n" \
  "\n" \
  PROG_INFO_OPTS_3 "\n" \
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
  "  sep/2008 Created by adaptation of {pnmfftfilter.c}.\n" \
  "  dec/2012 Made \"-zeroMean\" apply to the text output too.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmspectrum_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

/* We need to set these in order to get {isnan}. What a crock... */
#undef __STRICT_ANSI__
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <fftw3.h>

#include <jsfile.h>
#include <vec.h>
#include <r2.h>
#include <uint16_image.h>
#include <float_image.h>
#include <sample_conv.h>
#include <float_image_hartley.h>
#include <spectrum_table_exact.h>
#include <spectrum_table_binned.h>
#include <spectrum_table_convert.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <jspnm.h>
#include <argparser.h>

typedef struct options_t
  { /* Input image attributes: */
    string_vec_t inputImage;  /* Input image file name, or "-". */
    /* General parameters: */
    bool_t zeroMean;          /* TRUE to set the mean value to zero in the power spectrum. */
    bool_t averageImages;     /* TRUE to average all input images before computing the spectrum. */
    bool_t vignette;          /* TRUE to apply a vignetting mask before computing the spectrum. */
    bool_t symmetric;         /* True to output a symmetric spectrum. */
    /* Output image attributes: */
    char *outputImage;        /* output image file name, or "-", or NULL if none. */
    bool_t center;            /* TRUE to center the image on freq {(0,0)}. */
    float maxPower;           /* Nominal max power for sample conversion. */
    float minPower;           /* Nominal min power for logscale conversion, or 0.0 if linscale. */
    uint16_t maxPixel;        /* Output maxval requested by user, or 0 if not given. */
    /* output table attributes: */
    char *outputTable;        /* Name of output table file, or "-", or NULL if none. */
    int32_t freqRanges;           /* Number of ranges for binned, or 0 for exact. */
    /* Debugging options: */
    bool_t verbose;           /* TRUE to print diagnostic messages. */
  } options_t;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);

options_t *get_options(int32_t argc, char **argv);

float_image_t *read_image
  ( FILE *rd,
    int32_t *colsP,
    int32_t *rowsP,
    int32_t *chnsP,
    uint16_t *maxvalP,
    bool_t verbose
  );
  /* Reads a PBM/PGM/PPM image file from {rd}, converts it to a float
    image with samples in the range [0_1]. Returns the relevant image
    data. If {verbose} is true, prints image statistics to
    {stderr}. */

void accum_image(float_image_t *img_in, float_image_t *img_ot);
  /* Adds {img_in} to {img_ot}.  The two images  must have the same size on all alxes. */

void convert_to_spectrum(float_image_t *img_io, bool_t vignette);
  /* Replaces {img_io} by its (Hartley) power spectrum. If {vignette} is true, applies
    a vignetting window first. */

void apply_vignette(float_image_t *img_io, bool_t hann);
  /* Multiplies {img_io} by a Gaussian or Hann vignetting window. */

void write_image
  ( FILE *wr,            /* Where to write the output image. */
    float_image_t *fim,  /* Power spectrum image. */
    float v_min,         /* Nominal min power for logscale, 0 for linscale. */
    float v_max,         /* Nominal max power for logscale or linscale. */
    uint16_t maxval,     /* Max sample value in output PNM image. */
    bool_t verbose
  );
  /* Writes the float image {fim} to {wr} as a PBM/PGM/PPM image file.
    The output image file will have samples in {0..maxval}. 
    
    If {v_min} is 0, uses linear scale; if {v_min} is positive,
    uses logscale conversion with {v_min} as the reference value.

    The parameter {v_max} is the is image sample value that is to be mapped to {maxval};
    if it is 0 or negative, uses the actual max sample value in the image.
    
    If {verbose}
    is true, prints image statistics to {stderr}. */

void write_spectrum_table(FILE *wr, float_image_t *fim, int32_t nRanges, bool_t verbose);
  /* Writes to {wr} the spectrum power table, given the Hartley power spectrum
    image {fim}. If {nRanges} is zero, writes an exact table.
    If {nRanges} is positive, writes a binned table with {nRanges} bins.
    If {verbose}, prints diagnostics to {stderr}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    /* Parse the command line options: */
    options_t *o = get_options(argc, argv);

    /* Choose the output maxval: */
    uint16_t maxval_ot = (o->maxPixel > 0 ? o->maxPixel : 65535);
    
    /* Imput image dimensions: */
    int32_t chns = -1;
    int32_t cols = -1;
    int32_t rows = -1;
    
    /* Output accumulated power spectrum image (if any): */
    float_image_t *img_ot = NULL; /* Allocated after size is known. */
      
    int32_t NI = o->inputImage.ne; /* Number of images. */
    
    int32_t i;
    for (i = 0; i < NI; i++)
      { 
        /* Read the input image, get dimensions: */
        uint16_t maxval_in;
        FILE *rd = open_read(o->inputImage.e[i], o->verbose);
        float_image_t *img_in = read_image(rd, &cols, &rows, &chns, &maxval_in, o->verbose);

        /* Allocate output image if needed: */
        if (img_ot == NULL)
          { img_ot = float_image_new(chns, cols, rows);
            float_image_fill(img_ot, 0.0);
          }
          
        /* Accumulate power spectrum on output image: */
        if (o->averageImages)
          { accum_image(img_in, img_ot); }
        else
          { convert_to_spectrum(img_in, o->vignette);
            accum_image(img_in, img_ot);
          }
        
        float_image_free(img_in);
      }

    /* Convert total image to average: */
    if (NI != 1)
      { int32_t c;
        for (c = 0; c < chns; c++)
          { float_image_rescale_samples(img_ot, c, 0.0, (float)NI, 0.0, 1.0); }
      }

    if (o->averageImages)
      { /* Convert the input image average to the spectrum image: */
        convert_to_spectrum(img_ot, o->vignette);
      }
      
    if (o->zeroMean) { float_image_fill_pixel(img_ot, 0, 0, 0.0); }

    if (o->outputTable != NULL)
      { /* Write the filtered image, scaled and quantized: */
        FILE *wr = open_write(o->outputTable, TRUE);
        write_spectrum_table(wr, img_ot, o->freqRanges, o->verbose);
      }

    if (o->outputImage != NULL)
      { /* Write the power spectrum image: */
        if (o->symmetric)
          { /* Replace entries with same freq, opposite phase by their average: */
            int32_t c, x0, y0;
            for (c = 0; c < chns; c++)
              { for (x0 = 0; x0 <= cols/2; x0++)
                  { int32_t x1 = (cols - x0) % cols;
                    assert(x0 <= x1);
                    for (y0 = 0; y0 < rows; y0++)
                      { int32_t y1 = (rows - y0) % rows;
                        if ((x0 < x1) || ((x0 == x1) && (y0 < y1)))
                          { float *p0 = float_image_get_sample_address(img_ot, c, x0, y0);
                            float *p1 = float_image_get_sample_address(img_ot, c, x1, y1);
                            double v = (((double)(*p0)) + ((double)(*p1)))/2;
                            (*p0) = (float)v;
                            (*p1) = (float)v;
                          }
                      }
                  }
              }
          }
        if (o->center)
          { /* Shift the image so that {(0,0)} goes to the center: */
            int32_t c;
            for (c = 0; c < chns; c++) { float_image_shift(img_ot, c, cols/2, rows/2); }
          }
        FILE *wr = open_write(o->outputImage, TRUE);
        write_image(wr, img_ot, o->minPower, o->maxPower, maxval_ot, o->verbose);
        fclose(wr);
      }

    float_image_free(img_ot);
    if (o->verbose) { fprintf(stderr, "done.\n"); }
    return 0;
  }

float_image_t *read_image
  ( FILE *rd,
    int32_t *colsP,
    int32_t *rowsP,
    int32_t *chnsP,
    uint16_t *maxvalP,
    bool_t verbose
  )
  { uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    bool_t isMask = FALSE;
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, TRUE, verbose);
    if ((*colsP) < 0) { (*colsP) = pim->cols; } else { demand(pim->cols == (*colsP), "inconsistent image cols"); }
    if ((*rowsP) < 0) { (*rowsP) = pim->rows; } else { demand(pim->rows == (*rowsP), "inconsistent image rows"); }
    if ((*chnsP) < 0) { (*chnsP) = pim->chns; } else { demand(pim->chns == (*chnsP), "inconsistent image chns"); }
    (*maxvalP) = pim->maxval;
    uint16_image_free(pim);
    return fim;
  }

void accum_image(float_image_t *img_in, float_image_t *img_ot)
  {
    /* Get image dimensions: */
    int32_t NC = (int32_t)img_in->sz[0];
    int32_t NX = (int32_t)img_in->sz[1];
    int32_t NY = (int32_t)img_in->sz[2];

    /* Accumulate power spectrum on output image: */
    assert(img_ot->sz[0] == NC);
    assert(img_ot->sz[1] == NX);
    assert(img_ot->sz[2] == NY);
    int32_t c, x, y;
    for (c = 0; c < NC; c++) 
      { for (x = 0; x < NX; x++) 
          { for (y = 0; y < NY; y++)
              { double pv = float_image_get_sample(img_in, c, x, y);
                float *ov = float_image_get_sample_address(img_ot, c, x, y);
                (*ov) += (float)pv;
              }
          }
      }
  }
  
void convert_to_spectrum(float_image_t *img_io, bool_t vignette)
  {
    /* Get image dimensions: */
    int32_t NC = (int32_t)img_io->sz[0];
    int32_t NX = (int32_t)img_io->sz[1];
    int32_t NY = (int32_t)img_io->sz[2];

    /* Allocate the temp spectrum image: */
    float_image_t *img_pw = float_image_new(NC, NX, NY);

    if (vignette) 
      { bool_t hann = TRUE;
        apply_vignette(img_io, hann);
      }
    
    /* Compute the Hartley transform: */
    float_image_hartley_transform(img_io, img_pw);

    /* Store power spectrum on given image: */
    int32_t c, x, y;
    for (c = 0; c < NC; c++) 
      { for (x = 0; x < NX; x++) 
          { for (y = 0; y < NY; y++)
              { double pv = float_image_get_sample(img_pw, c, x, y);
                float *ov = float_image_get_sample_address(img_io, c, x, y);
                (*ov) = (float)(pv*pv);
              }
          }
      }
      
    float_image_free(img_pw);
  }
  
void apply_vignette(float_image_t *img_io, bool_t hann)
  {
    /* Get image dimensions: */
    int32_t NC = (int32_t)img_io->sz[0];
    int32_t NX = (int32_t)img_io->sz[1];
    int32_t NY = (int32_t)img_io->sz[2];
    
    double sigma = 0.5/3.0;
    
    auto double gauss_win(int32_t z, int32_t N);
    /* Gaussian mask with deviation {sigma*N}. */
    
    auto double hann_win(int32_t z, int32_t N);
    /* Hann mask for image size {N}. */
    
    double gauss_win(int32_t z, int32_t N)
      { double r = ((z + 0.5)/N - 0.5)/sigma;
        double w = exp(-0.5*r*r);
        return w;
      }

    double hann_win(int32_t z, int32_t N)
      { double r = (z + 0.5)/(double)N;
        double w = 0.5*(1.0 - cos(2*M_PI*r));
        return w;
      }

    /* Multiply each pixel by the mask: */
    int32_t c, x, y;
    for (x = 0; x < NX; x++) 
      { double wx = (hann ? hann_win(x, NX) : gauss_win(x, NX));
        for (y = 0; y < NY; y++)
          { double wy = (hann ? hann_win(y, NY) : gauss_win(y, NY));
            double w = wx*wy;
            for (c = 0; c < NC; c++) 
              { float *ov = float_image_get_sample_address(img_io, c, x, y);
                (*ov) *= (float)w;
              }
          }
      }
  }

void write_image
  ( FILE *wr,            /* Where to write the output image. */
    float_image_t *fim,  /* Power spectrum summed over all input images. */
    float v_min,         /* Nominal min power for logscale, 0 for linscale. */
    float v_max,         /* Nominal max power for logscale or linscale. */
    uint16_t maxval,     /* Max sample value in output PNM image. */
    bool_t verbose
  )
  { int32_t chns = (int32_t)fim->sz[0];
    int32_t c;
    demand(v_min >= 0, "invalid {v_min}");

    /* Ensure {v_max} is defined and positive: */
    if (v_max <= 0)
      { /* Find true {v_max} over all channels (but leave {v_min} unchanged): */
        float f_min = 0.0f, f_max = 1.0e-38f; /* To avoid division by zero if {fim} is all zeros. */
        for (c = 0; c < chns; c++)
          { float_image_update_sample_range(fim, c, &f_min, &f_max); }
        if (verbose) { fprintf(stderr, "max power = %13.5e\n", f_max); }
        v_max = f_max;
      }
    assert(v_max > 0.0);

    /* Convert image and {v_max} to logscale if so requested: */
    if (v_min > 0)
      { /* Apply log-scale conversion to the image: */
        if (verbose) { fprintf(stderr, "applying log scale [ %13.5e _  %13.5e ] -->", v_min, v_max); }
        double base = M_E;
        for (c = 0; c < chns; c++) 
          { float_image_log_scale(fim, c, v_min, base); }
        /* Apply log-scale conversion to {v_max} too: */
        double logBase = 1.0;
        v_max = fmaxf(1.0e-38f, sample_conv_log(v_max, v_min, logBase));
        /* The old {v_min} should be mapped to zero: */
        v_min = 0.0;
        if (verbose) { fprintf(stderr, " [ %13.5e _  %13.5e ]\n", v_min, v_max); }
      }

    /* At this point we must have {v_min == 0} in any case: */
    assert(v_min == 0.0);
    assert(v_max > 0.0);

    /* Write image with linear scaling {[0 _ v_max]} --> {[0 _ maxval}: */
    double vLo[chns];
    double vHi[chns];
    for (c = 0; c < chns; c++) { vLo[c] = 0; vHi[c] = v_max; }
    bool_t isMask = FALSE;
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, chns, vLo, vHi, NULL, maxval, TRUE, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    uint16_image_free(pim);
  }

void write_spectrum_table(FILE *wr, float_image_t *fim, int32_t nRanges, bool_t verbose)
  {
    int32_t chns = (int32_t)fim->sz[0];
    int32_t cols = (int32_t)fim->sz[1];
    int32_t rows = (int32_t)fim->sz[2];

    /* Create a binned spectrum table {tb}: */
    spectrum_table_binned_t tb;
    int32_t c;
    if (nRanges == 0)
      { /* Create an exact spectrum table {tx}: */
        if (verbose) { fprintf(stderr, "gathering the exact spectrum table...\n"); }
        spectrum_table_exact_t tx = spectrum_table_exact_new(0);
        if (verbose) { fprintf(stderr, "  %12u raw entries\n", tx.ne); }
        for (c = 0; c < chns; c++)
          { spectrum_table_exact_append_all(fim, c, &tx, FALSE); }
        spectrum_table_exact_sort(&tx, verbose);
        if (verbose) { fprintf(stderr, "  %12u sorted entries\n", tx.ne); }
        /* Convert {tx} to a binned table {tb}, preserving the exactness as fas as possible: */
        tb = spectrum_table_convert_exact_to_binned(&tx, cols, rows);
        if (verbose) { fprintf(stderr, "  %12u binned entries\n", tb.ne); }
        spectrum_table_exact_trim(&tx, 0);
      }
    else
      { /* Create a binned spectrum table, put it in {tb}: */
        if (verbose) { fprintf(stderr, "gathering the binned spectrum table...\n"); }
        tb = spectrum_table_binned_make(nRanges);
        if (verbose) { fprintf(stderr, "  %12u binned entries\n", tb.ne); }
        for (c = 0; c < chns; c++)
          { spectrum_table_binned_add_all(fim, c, &tb, FALSE); }
      }

    /* Write the binned spectrum table to the output: */
    int32_t k;
    double fprev = 0.0;
    char *tbfmt = "%14.8f %14.10f %14.8f  %10.0f  %18.12f\n";
    for (k = 0; k < tb.ne; k++)
      { spectrum_table_binned_entry_t *ek = &(tb.e[k]);
        assert(fprev == ek->fmin);
        assert(ek->fmin <= ek->fmid);
        assert(ek->fmid <= ek->fmax);
        fprintf(wr, tbfmt, ek->fmin, ek->fmid, ek->fmax, ek->nTerms, ek->power);
        fprev = ek->fmax;
      }
    fflush(wr);
    spectrum_table_binned_trim(&tb, 0);
  }

#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

#define MAX_SIZE (32*1024)
  /* A limit on image size, to avoid humongous mallocs. */

#define MAX_FREQ_RANGES (1024*1024)
  /* A limit on spectrun table size, to avoid humongous mallocs. */

options_t *get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);

    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");

    /* Arguments related to image output: */

    o->zeroMean = argparser_keyword_present(pp, "-zeroMean");

    o->averageImages = argparser_keyword_present(pp, "-averageImages");

    o->vignette = argparser_keyword_present(pp, "-vignette");

    if (argparser_keyword_present(pp, "-outputImage"))
      { if (argparser_keyword_present_next(pp, "NONE"))
          { o->outputImage = NULL; }
        else
          { o->outputImage = argparser_get_next(pp); }
      }
    else
      { o->outputImage = NULL; }

    o->center = argparser_keyword_present(pp, "-center");

    o->symmetric = argparser_keyword_present(pp, "-symmetric");

    if (argparser_keyword_present(pp, "-scale"))
      { if (argparser_keyword_present_next(pp, "linear"))
          { o->minPower = 0.0; /* Implies `use linear scale'. */ }
        else if (argparser_keyword_present_next(pp, "log"))
          { o->minPower = (float)argparser_get_next_double(pp, 1.0e-38, 1.0e+38); }
        else
          { argparser_error(pp, "invalid scale type"); }
      }
    else
      { o->minPower = 0.0; /* Implies `use linear scale'. */ }

    if (argparser_keyword_present(pp, "-maxPower"))
      { if (argparser_keyword_present_next(pp, "AUTO"))
          { o->maxPower = 0.0; /* Implies `use actual max power'. */ }
        else
          { o->maxPower = (float)argparser_get_next_double(pp, 1.0e-38, 1.0e+38); }
      }
    else
      { o->maxPower = 0.0; /* Implies `use actual max power'. */ }

    if (argparser_keyword_present(pp, "-maxPixel"))
      { o->maxPixel = (uint16_t)argparser_get_next_int(pp, 1, PNM_MAX_SAMPLE); }
    else
      { o->maxPixel = 0; /* Implies `use {max(255,maxval_in)}'. */ }

    /* Arguments related to tabular output: */

    if (argparser_keyword_present(pp, "-outputTable"))
      { if (argparser_keyword_present_next(pp, "NONE"))
          { o->outputTable = NULL; }
        else
          { o->outputTable = argparser_get_next(pp); }
      }
    else
      { o->outputTable = NULL; }

    if (argparser_keyword_present(pp, "-freqRanges"))
      { o->freqRanges = (int32_t)argparser_get_next_int(pp, 1, MAX_FREQ_RANGES); }
    else if (argparser_keyword_present(pp, "-exact"))
      { o->freqRanges = 0; /* Implies `exact output'. */ }
    else 
      { o->freqRanges = 0; /* Implies `exact output'. */ }

    if ((o->outputTable == NULL) && (o->outputImage == NULL))
      { argparser_error(pp, "Must use either \"-outputImage\" or \"-outputTable\", or both."); }

    /* Other keyword arguments: */

    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Go to positional arguments: */
    argparser_skip_parsed(pp);

    /* Parse optional input file name: */
    o->inputImage = string_vec_new(10);
    int32_t NI = 0;
    if (argparser_next(pp) == NULL)
      { 
        o->inputImage.e[0] = "-";
        NI = 1;
      }
    else
      { 
        while (argparser_next(pp) != NULL)
          { string_vec_expand(&(o->inputImage), NI);
            o->inputImage.e[NI] = argparser_get_next(pp);
            NI++;
          }
      }
    string_vec_trim(&(o->inputImage), NI);

    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }

