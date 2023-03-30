#define PROG_NAME "fni_to_pnm"
#define PROG_DESC "convert a float-valued FNI image file to a PGM or PPM file"
#define PROG_VERS "1.0"

/* Last edited on 2023-03-29 12:27:42 by stolfi */

#define PROG_C_COPYRIGHT "Copyright © 2005 by the State University of Campinas (UNICAMP)."

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -channel {CHAN} | -channels {CHANR} {CHANG} {CHANB} ] \\\n" \
  "    [ " pst_scaling_parse_min_any_HELP " ] \\\n" \
  "    [ " pst_scaling_parse_max_any_HELP " ] \\\n" \
  "    [ " pst_scaling_parse_center_any_HELP " ] \\\n" \
  "    [ " pst_scaling_parse_width_any_HELP " ] \\\n" \
  "    [ " pst_scaling_parse_uniform_HELP " ] \\\n" \
  "    [ -maxval {MAXVAL} ] \\\n" \
  "    [ -isMask {ISMASK} ] \\\n" \
  "    [ -nanval {NANVAL} ] \\\n" \
  "    < {FNI_FILE} \\\n" \
  "    > {PNM_FILE}"

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
  "  pnm_to_fni(1), {float_image.h}.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2005-12-10 by Jorge Stolfi, UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  aug/2009 Replacement of NaNs added by R. F. V. Saracchini.\n" \
  "  aug/2009 \"-nanval\" option added by J. Stolfi.\n" \
  "  2010-08-14 \"-isMask\" flag added by J. Stolfi.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program reads a float-valued multichannel image" \
  " (FNI) file and writes either one of its channels" \
  " as a raw PGM file (magic number \"P5\"), or three of its" \
  " channels as a raw PPM file (magic number \"P6\").  See" \
  " {float_image.h} for a description of the FNI format.\n" \
  "\n" \
  "  The indices of the channels to be converted are specified" \
  " through the \"-channel\" option (for PGM output) or the" \
  " \"-channels\" option (for PPM output).  Channel indices start" \
  " at zero; channels that do not exist in the input" \
  " image are assumed to be all zeros.  If neither option is present," \
  " the input must have either one or three channels, and then" \
  " \"-channel 0\" or \"-channels 0 1 2\" is assumed, respectively.\n" \
  "\n" \
  "  In either case, the output PGM/PPM image has the same" \
  " dimensions as the input float image, and its samples" \
  " range from 0 to {MAXVAL} (default 65535).\n" \
  "\n" \
  "  Input sample values are clipped to a certain range" \
  " [{VMIN} _ {VMAX}], then mapped affinely (linearly)" \
  " from that range to {[0_1]}, with clipping; and then that range is encoded to" \
  " the integers {0..MAXVAL}, as specified by the \"-isMask\" option.  The" \
  " encoding is essentially linear, without any gamma mapping.  Any {NaN} samples" \
  " in the input file are mapped directly to the sample value defined" \
  " by the \"-nanval\" argument.  The top scanline of" \
  " the output image is row 0 of the input.\n" \
  "\n" \
  "SPECIFYING THE SAMPLE SCALING\n" \
  "  The input scaling range [{VMIN} _ {VMAX}] is specified by" \
  " the options " pst_scaling_option_list_INFO ".\n" \
  "\n" \
  "  " pst_scaling_use_actual_range_INFO "  " \
  pst_scaling_use_actual_range_with_uniform_INFO "\n" \
  "\n" \
  "  In any case, at most two of these four parameters may be specified. " \
  " " pst_scaling_complete_params_INFO ""
  
#define PROG_INFO_OPTS \
  "  -channel {CHAN}\n" \
  "    Requests a PGM output image containing channel {CHAN}" \
  " of the input image.\n" \
  "\n" \
  "  -channels {CHAN_R} {CHAN_G} {CHAN_B}\n" \
  "    Requests a PPM output image using channels {CHAN_R}," \
  " {CHAN_G}, and {CHAN_B} of the input image as the" \
  " Red, Green, and Blue components.  The three channel" \
  " indices need not be distinct.\n" \
  "\n" \
  "  " pst_scaling_parse_min_one_HELP "\n" \
  "  " pst_scaling_parse_min_RGB_HELP "\n" \
  "    " pst_scaling_parse_min_INFO \
  " of the input scaling range.  " \
  pst_double_vec_spec_den_INFO "  " \
  pst_scaling_num_values_INFO "\n" \
  "\n" \
  "  " pst_scaling_parse_max_one_HELP "\n" \
  "  " pst_scaling_parse_max_RGB_HELP "\n" \
  "    " pst_scaling_parse_max_INFO \
  " of the input scaling range.  " \
  pst_double_vec_spec_den_INFO "  " \
  pst_scaling_num_values_INFO "\n" \
  "\n" \
  "  " pst_scaling_parse_center_one_HELP "\n" \
  "  " pst_scaling_parse_center_RGB_HELP "\n" \
  "    " pst_scaling_parse_center_INFO \
  " of the input scaling range.  " \
  pst_double_vec_spec_den_INFO "  " \
  pst_scaling_num_values_INFO "\n" \
  "\n" \
  "  " pst_scaling_parse_width_one_HELP "\n" \
  "  " pst_scaling_parse_width_RGB_HELP "\n" \
  "    " pst_scaling_parse_width_INFO \
  " of the input scaling range.  " \
  pst_double_vec_spec_den_INFO "  " \
  pst_scaling_num_values_INFO "\n" \
  "\n" \
  pst_scaling_parse_uniform_HELP_INFO ".\n" \
  "\n" \
  "  -maxval {INTEGER}\n" \
  "    Specifies the maximum output sample value {MAXVAL}.\n" \
  "\n" \
  "  -isMask {ISMASK}\n" \
  "    This optional Boolean argument specifies the interpretation" \
  " of integer sample values in the output file, specifically" \
  " how the float values already normalized to the range {[0_1]} are mapped" \
  " to the integer values {0..MAXVAL}.  If {ISMASK} is true (\"T\" or 1)," \
  " " sample_conv_0_1_isMask_true_INFO "  If {ISMASK} is false (\"F\" or 0)," \
  " " sample_conv_0_1_isMask_false_INFO "  The default is \"-isMask F\".\n" \
  "\n" \
  "  -nanval {UNDVAL}\n" \
  "    This optional parameter specifies the floating-point value to" \
  " be substituted for {NAN} samples. After substitution, that value" \
  " will be converted by the same formula as the original samples.  If this" \
  " parameter is not specified, the program will use the midpoint of the" \
  " nominal input range of each channel, as specified by the scaling arguments."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <argparser.h>
#include <float_image.h>
#include <float_image_to_uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <sample_conv.h>
#include <bool.h>

#include <pst_scaling.h>
#include <pst_basic.h>

#define INF INFINITY

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { int NC;              /* Number of channels to output, or -1 if unknown. */
    uint16_t maxval; /* Maximum output sample value. */
    bool_t isMask;       /* Interpretation of integer sample values. */
    /* Input channel data (indexed {0..NC-1}): */
    int32_vec_t channel;   /* Indices of input channels to convert; empty vec if not given. */
    double_vec_t min;    /* Low endpoint of scaling range; empty vec if not given. */
    double_vec_t max;    /* High endpoint of scaling range; empty vec if not given. */
    double_vec_t ctr;    /* Center of scaling range; empty vec if not given. */
    double_vec_t wid;    /* Width of scaling range; empty vec if not given. */
    bool_t uniform;      /* Obtains default scaling args from whole image rather than single channel. */
    float nanval;        /* Value to be substituted for {NAN} samples, or {NAN} if none. */
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char** argv);

options_t *ftp_parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

float_image_t *ftp_read_fni_image(FILE *rd, int *NC, int *NX, int *NY);
  /* Reads a float-valued image from the FNI file {rd}.
    Sets {*NC,*NX,*NY} to the image dimensions. */

void ftp_float_image_write_pgm(FILE *wr, float_image_t *fim, bool_t isMask, int c, double lo, double hi, uint16_t maxval);
  /* Writes channel {c} of image {fim} to file {wr} as a PGM image,
    mapping each pixel from {[lo _ hi]} to {[0..maxval]}. Note that
    row 0 of {fim} is the *bottom* row of the PGM image. */

void ftp_write_ppm_image(FILE *wr, float_image_t *fim, bool_t isMask, int ch[], double lo[], double hi[], uint16_t maxval);
  /* Writes channels {ch[0..2]} of image {fim} to file {wr} as a PPM image,
    scaling each sample of channel {ch[k]} linearly from {[lo[k] _ hi[k]]} to
    {[0..maxval]}. Note that row 0 of {fim} is the *bottom* row of the PPM image. */

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = ftp_parse_options(argc, argv);

    /* Read image, get its depth: */
    int NCI, NXI, NYI; /* Input image dimensions. */
    float_image_t *fim = ftp_read_fni_image(stdin, &NCI, &NXI, &NYI);
    
    /* If the number of output channels was not specified by the user, use {NCI}. */
    if (o->NC < 0) { o->NC = NCI; }

    /* Make sure we have 1 or 3 channels: */
    if ((o->NC != 1) && (o->NC != 3))
      { /* The user omitted "-channel"/"-channels" for an image with weird depth: */
        fprintf(stderr, "input image has %d channels\n", o->NC);
        fprintf(stderr, "must specify \"-channel\" or \"-channels\"\n");
        exit(1);
      }

    /* Provide the default channel indices, if not specified: */
    pst_scaling_fix_channels(o->NC, &(o->channel));

    /* Complete the scaling parameters, with defaults and derived values: */
    pst_scaling_fix_params
      ( o->NC, o->uniform,
        &(o->min), &(o->max), &(o->ctr), &(o->wid),
        fim, &(o->channel)
      );
      
    /* Replace {NAN}s: */
    int c;
    for (c = 0; c < NCI; c++)
      { float v = o->nanval;
        if (isnan(v)) { v = (float)o->ctr.e[c]; }
        float_image_replace_nan_samples(fim, c, v);
      }

    /* Convert the image {fim} to PGM or PPM and write the result: */
    if (o->NC == 1)
      { ftp_float_image_write_pgm(stdout, fim, o->isMask, o->channel.e[0], o->min.e[0], o->max.e[0], o->maxval); }
    else if (o->NC == 3)
      { ftp_write_ppm_image(stdout, fim, o->isMask, o->channel.e, o->min.e, o->max.e, o->maxval); }
    else
      { assert(FALSE); }

    return 0;
  }

float_image_t *ftp_read_fni_image(FILE *rd, int *NC, int *NX, int *NY)
  { float_image_t *fim = float_image_read(rd);
    (*NC) = (int)fim->sz[0];
    (*NX) = (int)fim->sz[1];
    (*NY) = (int)fim->sz[2];
    return fim;
  }

void ftp_float_image_write_pgm(FILE *wr, float_image_t *fim, bool_t isMask, int c, double lo, double hi, uint16_t maxval)
  { int ch[1]; /* Channels to convert. */
    ch[0] = c;
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, 1, &lo, &hi, ch, maxval, TRUE, TRUE);
    /* Write to disk: */
    bool_t forceplain = FALSE;
    bool_t verbose = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    fflush(wr);
    /* Cleanup: */
    uint16_image_free(pim);
  }

void ftp_write_ppm_image(FILE *wr, float_image_t *fim, bool_t isMask, int ch[], double lo[], double hi[], uint16_t maxval)
  { /* Convert {fim} to PPM image: */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, 3, lo, hi, ch, maxval, TRUE, TRUE);
    /* Write to disk: */
    bool_t forceplain = FALSE;
    bool_t verbose = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    fflush(wr);
    /* Cleanup: */
    uint16_image_free(pim);
  }

options_t *ftp_parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)malloc(sizeof(options_t));

    /* Get number of channels and channel list (if specified): */
    o->NC = -1;  /* Unknown number of channels. */
    if (argparser_keyword_present(pp, "-channel"))
      { o->NC = 1; }
    else if (argparser_keyword_present(pp, "-channels"))
      { o->NC = 3; }
    if (o->NC >= 0)
      { o->channel = pst_int32_vec_parse(pp, &(o->NC)); }
    else
      { o->channel = int32_vec_new(0); }

    /* Parse input range specs: */
    o->min = pst_scaling_parse_range_option(pp, "-min",    &(o->NC));
    o->max = pst_scaling_parse_range_option(pp, "-max",    &(o->NC));
    o->ctr = pst_scaling_parse_range_option(pp, "-center", &(o->NC));
    o->wid = pst_scaling_parse_range_option(pp, "-width",  &(o->NC));

    /* Uniform scaling option: */
    o->uniform = pst_scaling_parse_uniform(pp, FALSE);

    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (uint16_t)argparser_get_next_int(pp, 1, 65535); }
    else
      { o->maxval = 65535; }

    if (argparser_keyword_present(pp, "-isMask"))
      { o->isMask = argparser_get_next_bool(pp); }
    else
      { o->isMask = FALSE; }

    if (argparser_keyword_present(pp, "-nanval"))
      { o->nanval = (float)argparser_get_next_double(pp, -INF, +INF); }
    else
      { o->nanval = NAN; }

    argparser_finish(pp);

    return o;
  }
