#define PROG_NAME "fni_to_pnm"
#define PROG_DESC "convert a float-valued FNI image file to a PGM or PPM file"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-22 18:48:41 by stolfi */

#define PROG_C_COPYRIGHT "Copyright © 2005 State University of Campinas (UNICAMP).  Run \"" PROG_NAME " -info\" for details"

/* !!! Add "-gamma" option !!! */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -channel {CHAN} | -channels {CHANR} {CHANG} {CHANB} ] \\\n" \
  "    " sample_scaling_options_parse_HELP " \\\n" \
  "    [ -maxval {MAXVAL} ] \\\n" \
  "    [ -isMask {ISMASK} ] \\\n" \
  "    [ -nanval {NANVAL} ] \\\n" \
  "    " imgc_parse_y_axis_HELP " \\\n" \
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
  "  2024-12-24 Added \"-yAxis\". J.Stolfi.\n" \
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
  " by the \"-nanval\" argument.\n" \
  "\n" \
  "  The \"-yAxis\" option determines whether row 0 of the" \
  " FNI image is the top or bottom row of the PNM image.\n" \
  "\n" \
  "SPECIFYING THE SAMPLE SCALING\n" \
  "  " sample_scaling_options_parse_HELP_INFO ""
  
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
  sample_scaling_options_parse_HELP_INFO ".\n" \
  "\n" \
  imgc_parse_y_axis_INFO_OPTS(imgc_parse_y_axis_INFO_OPTS_default_pbm) "\n" \
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
  "  -nanval {NANVAL}\n" \
  "    This optional parameter specifies the floating-point value to" \
  " be substituted for {NAN} samples. After substitution, that value" \
  " will be converted by the same formula as the original samples.  If this" \
  " parameter is not specified, the program will use the midpoint of the" \
  " nominal input range of each channel, as specified by the scaling arguments."

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <argparser.h>
#include <argparser_extra.h>
#include <float_image.h>
#include <float_image_to_uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <sample_conv.h>
#include <image_coords.h>
#include <bool.h>
#include <sample_scaling.h>

#include <pst_basic.h>

#define INF INFINITY

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { int32_t NC;            /* Number of channels to output, or -1 if unknown. */
    uint16_t maxval;       /* Maximum output sample value. */
    bool_t isMask;         /* Interpretation of integer sample values. */
    /* Input channel data (indexed {0..NC-1}): */
    int32_vec_t channel;   /* Indices of input channels to convert; empty vec if not given. */
    sample_scaling_options_t scaling; 
    float nanval;        /* Value to be substituted for {NAN} samples, or {NAN} if none. */
    bool_t yUp;          /* If true, row 0 of the FNI image is bottom of PNM. */
  } options_t;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char** argv);

options_t *ftp_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

float_image_t *ftp_read_fni_image(FILE *rd, int32_t *NC, int32_t *NX, int32_t *NY);
  /* Reads a float-valued image from the FNI file {rd}.
    Sets {*NC,*NX,*NY} to the image dimensions. */

void ftp_float_image_write_pgm(FILE *wr, float_image_t *fim, bool_t isMask, int32_t c, double lo, double hi, uint16_t maxval, bool_t yUp);
  /* Writes channel {c} of image {fim} to file {wr} as a PGM image,
    mapping each pixel from {[lo _ hi]} to {[0..maxval]}. Note that
    row 0 of {fim} is the *bottom* row of the PGM image. */

void ftp_write_ppm_image
  ( FILE *wr,
    float_image_t *fim,
    bool_t isMask,
    int32_t ch[],
    double lo[],
    double hi[],
    uint16_t maxval,
    bool_t yUp
  );
  /* Writes channels {ch[0..2]} of image {fim} to file {wr} as a PPM image,
    scaling each sample of channel {ch[k]} linearly from {[lo[k] _ hi[k]]} to
    {[0..maxval]}. Note that row 0 of {fim} is the *bottom* row of the PPM image. */
  /* !!! What happens to {NAN}? !!! */
    
/* !!! Add {nanval} option !!! */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = ftp_parse_options(argc, argv);

    /* Read image, get its depth: */
    int32_t NCI, NXI, NYI; /* Input image dimensions. */
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
    sample_scaling_fix_channels(o->NC, &(o->channel));

    /* Complete the scaling parameters, with defaults and derived values: */
    sample_scaling_fix_params(&(o->scaling), &(o->channel), o->NC, fim);
      
    /* Replace {NAN}s by mid-range sample values: */
    for (int32_t c = 0; c < NCI; c++)
      { float v = o->nanval;
        if (isnan(v)) { v = (float)o->scaling.ctr.e[c]; }
        float_image_replace_nan_samples(fim, c, v);
      }

    /* Convert the image {fim} to PGM or PPM and write the result: */
    if (o->NC == 1)
      { int32_t ch = o->channel.e[0];
        double vmin = o->scaling.min.e[0];
        double vmax = o->scaling.max.e[0];
        ftp_float_image_write_pgm(stdout, fim, o->isMask, ch, vmin, vmax, o->maxval, o->yUp); }
    else if (o->NC == 3)
      { 
        int32_t *chans = o->channel.e;
        double *vmin = o->scaling.min.e;
        double *vmax = o->scaling.max.e;
        ftp_write_ppm_image(stdout, fim, o->isMask, chans, vmin, vmax, o->maxval, o->yUp); }
    else
      { assert(FALSE); }

    return 0;
  }

float_image_t *ftp_read_fni_image(FILE *rd, int32_t *NC, int32_t *NX, int32_t *NY)
  { float_image_t *fim = float_image_read(rd);
    (*NC) = (int32_t)fim->sz[0];
    (*NX) = (int32_t)fim->sz[1];
    (*NY) = (int32_t)fim->sz[2];
    return fim;
  }

void ftp_float_image_write_pgm
  ( FILE *wr,
    float_image_t *fim,
    bool_t isMask,
    int32_t c,
    double lo,
    double hi,
    uint16_t maxval,
    bool_t yUp
  )
  { int32_t ch[1]; /* Channels to convert. */
    ch[0] = c;
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, 1, &lo, &hi, ch, maxval, yUp, TRUE);
    /* Write to disk: */
    bool_t forceplain = FALSE;
    bool_t verbose = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    fflush(wr);
    /* Cleanup: */
    uint16_image_free(pim);
  }

void ftp_write_ppm_image
  ( FILE *wr,
    float_image_t *fim,
    bool_t isMask,
    int32_t ch[],
    double lo[],
    double hi[],
    uint16_t maxval,
    bool_t yUp
  )
  { /* Convert {fim} to PPM image: */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, 3, lo, hi, ch, maxval, yUp, TRUE);
    /* Write to disk: */
    bool_t forceplain = FALSE;
    bool_t verbose = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    fflush(wr);
    /* Cleanup: */
    uint16_image_free(pim);
  }

options_t *ftp_parse_options(int32_t argc, char **argv)
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
      { o->channel = argparser_get_next_int32_vec(pp, &(o->NC)); }
    else
      { o->channel = int32_vec_new(0); }

    /* Parse input range specs: */
    o->scaling = sample_scaling_parse_options(pp, &(o->NC));

    o->yUp = FALSE;
    imgc_parse_y_axis(pp, &(o->yUp));

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
