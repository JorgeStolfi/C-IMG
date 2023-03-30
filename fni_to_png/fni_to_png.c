#define PROG_NAME "fni_to_pnm"
#define PROG_DESC "convert a float-valued FNI image file to a PNG file"
#define PROG_VERS "1.0"

/* Last edited on 2023-03-29 12:51:30 by stolfi */

#define PROG_C_COPYRIGHT "Copyright © 2021 by the State University of Campinas (UNICAMP)."

/* !!! Add "-gamma" option !!! */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -channels {CHAN}... ] \\\n" \
  "    [ -min {VMIN}... ] \\\n" \
  "    [ -max {VMAX}... ] \\\n" \
  "    [ -maxval {MAXVAL} ] \\\n" \
  "    [ -isMask {ISMASK} ] \\\n" \
  "    [ -nanval {NANVAL} ] \\\n" \
  "    [ -verbose {VERB} ] \\\n" \
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
  "  Created 2021-08-24 by Jorge Stolfi, UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2021-08-24 Created based on {fni_to_pnm.c}. J. Stolfi.\n" \
  "  2021-08-24 Added \"-verbose\" option. J. Stolfi.\n" \
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
  " (FNI) file and writes a subset of one to four of its channels" \
  " as a PNG file.  See {float_image.h} for a description of the FNI format.\n" \
  "\n" \
  "  The indices of the channels to be converted are specified" \
  " through the \"-channels\" option.The type of the output image depends on the number of selected channels \n" \
  "\n" \
  "    1 channel: grayscale.\n" \
  "    2 channels: grayscale + opacity.\n" \
  "    3 channels: RGB.\n" \
  "    3 channels: RGB + opacity.\n" \
  "\n" \
  "  In any case, the output PNG image has the same" \
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
  " the output image is row 0 of the input."
  
#define PROG_INFO_OPTS \
  "  -channel {CHAN}\n" \
  "  -channels {CHAN}...\n" \
  "    These two mutually exclusive optional arguments pecify the" \
  " indices of between 1 and 4 channels to use in the output image.  Channel" \
  " indices start at zero; channels that do not exist in the input" \
  " image are assumed to be all zeros.  If neither option is present," \
  " the input must have at most 4 channels.  The channel" \
  " indices need not be distinct.\n" \
  "\n" \
  "  -min {VMIN}...\n" \
  "  -max {VMAX}...\n" \
  "    These optional arguments specify the range of values from the input image that" \
  " is to be mapped to {0..MAXVAL} in the output.  Each keyword should be followed by" \
  " one real value, that applies to all channels, or by one value for each selected channel.  If" \
  " either parameter is omitted, or the value given is \"nan\", the program uses the minimum or maxmium" \
  " sample value in the respective channel.\n" \
  "\n" \
  "  -uniform {UNI_FLAG}.\n" \
  "    This optional parameter is relevant only if {VMIN} or {VMAX} is unspecified for some" \
  " channel(s).  If it is present and true, the maxmimum and mnimum sample" \
  " values are computed for all unspecified channels, instead of separately for each" \
  " channel (the default).\n" \
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
  " nominal input range of each channel, as specified by the scaling arguments." \
  "\n" \
  "  -verbose {VERB}\n" \
  "    This optional parameter requests the printout of image statistics and" \
  " conversion parameters.  The {VERB} value may be \"T\" or 1 to get the" \
  " diagnostics, or \"F\" of 0 to suppress them.  If omitted, the program" \
  " assumes \"-verbose F\"."

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <argparser.h>
#include <float_image.h>
#include <float_image_to_uint16_image.h>
#include <uint16_image_write_png.h>
#include <uint16_image.h>
#include <sample_conv.h>
#include <bool.h>

#include <pst_scaling.h>
#include <pst_basic.h>

#define INF INFINITY

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { uint16_t maxval;         /* Maximum output sample value. */
    bool_t isMask;           /* Interpretation of integer sample values. */
    int32_vec_t channels;      /* Indices of input channels to convert; empty vec if not given. */
    double_vec_t min;        /* Low endpoint of scaling range; empty vec if not given. */
    double_vec_t max;        /* High endpoint of scaling range; empty vec if not given. */
    bool_t uniform;          /* Obtains default scaling args from all unspec channels rather than each channel. */
    float nanval;            /* Value to be substituted for {NAN} samples, or {NAN} if none. */
    bool_t verbose;          /* True to print image stats and conversion info. */
  } options_t;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char** argv);

options_t *ftp_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

double_vec_t ftp_parse_double_vec(argparser_t *pp, char *key);
int32_vec_t ftp_parse_int32_vec(argparser_t *pp, char *key);
  /* If keyword {key} is present, parses one or more numbers after it. */

void ftp_fix_scaling_param_vec(int32_t nc, double_vec_t *vex);
  /* If {vex} has {nc} elements, does nothing. If {vex} is an empty vector, provides a vector of {nc} {NAN}s.
    If {vex} has one value and {nc>1}, replicates that value to all other elements. Fails othewise. */
    
void ftp_fill_missing_range_params(float_image_t *fim, int32_t nc, int32_t ch[], bool_t uniform, double_vec_t *vmin, double_vec_t *vmax);
  /* If any of the elements of {vmin} and/or {vmax} is {NAN}, supplies a value from the actual range
    of sample in the corresponding channel(s) of {fim}.  Ignores elements {vmin.e[c]} or {vmax.e[c]}
    such that {ch[c]} is not a valid {fim} channel index. */
    
float_image_t *ftp_read_fni_image(FILE *rd, int32_t *NC, int32_t *NX, int32_t *NY);
  /* Reads a float-valued image from the FNI file {rd}.
    Sets {*NC,*NX,*NY} to the image dimensions. */

void ftp_write_png_image
  ( FILE *wr, 
    float_image_t *fim, 
    bool_t isMask, 
    int32_t nc, 
    int32_t ch[], 
    double lo[], 
    double hi[],
     uint16_t maxval, 
    float nanval, 
    bool_t verbose
  );
  /* Writes channels {ch[0..nc-1]} of image {fim} to file {wr} as a PNG image,
    scaling each sample of channel {ch[k]} linearly from {[lo[k] _ hi[k]]} to
    {[0..maxval]}.  Note that row 0 of {fim} is the *bottom* row of the PNG image.
    
    Any {NAN} samples are replaced by {nanval} before conversion. If {nanval} itself
    is {NAN}, it defaults to the midpoint of the channel's range {[lo[c]__hi[c]]}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t *o = ftp_parse_options(argc, argv);

    /* Read image, get its depth: */
    int32_t NCI, NXI, NYI; /* Input image dimensions. */
    float_image_t *fim = ftp_read_fni_image(stdin, &NCI, &NXI, &NYI);
    
    /* If the number of output channels was not specified by the user, use {NCI}. */
    if (o->channels.ne == 0)
      { int32_vec_expand(&(o->channels), NCI-1);
        for (int32_t c = 0; c < NCI; c++) { o->channels.e[c] = c; }
      }
    int32_t nc = o->channels.ne;
    int32_t *ch = o->channels.e;
    if ((nc < 1) || (nc > 4))
      { fprintf(stderr, "** number of channels %d should be in 1..4\n", nc);
        demand(FALSE, "aborted");
      }

    /* Complete the scaling parameters, with defaults and derived values: */
    ftp_fix_scaling_param_vec(nc, &(o->min));
    ftp_fix_scaling_param_vec(nc, &(o->max));
    ftp_fill_missing_range_params(fim, nc, ch, o->uniform, &(o->min), &(o->max));
    
    /* Convert the image {fim} to PNG and write the result: */
    ftp_write_png_image(stdout, fim, o->isMask, nc, ch, o->min.e, o->max.e, o->maxval, o->nanval, o->verbose);

    return 0;
  }

void ftp_fix_scaling_param_vec(int32_t nc, double_vec_t *vex)
  { assert((nc >= 1) && (nc <= 4));
    if (vex->ne == nc)
      { /* Nothing to do: */
      }
    else if (vex->ne == 0)
      { /* Set to a vecor of {NAN}: */
        (*vex) = double_vec_new(nc);
        for (int32_t c = 0; c < nc; c++) { vex->e[c] = NAN; }
      }
    else if ((vex->ne == 1) && (nc > 1))
      { /* Replicate to all elems: */
        double_vec_expand(vex, nc-1);
        for (int32_t c = 1; c < nc; c++) { vex->e[c] = vex->e[0]; }
      }
    else
      { demand(FALSE, "invalid range spec"); }
  }
      
void ftp_fill_missing_range_params(float_image_t *fim, int32_t nc, int32_t ch[], bool_t uniform, double_vec_t *vmin, double_vec_t *vmax)
  { int32_t NCI = (int32_t)fim->sz[0];
    float smin, smax;
    if (uniform)
      { /* Find global range of all channels with unspecified range: */
        smin = +INF; smax = -INF;
        for (int32_t c = 0; c < nc; c++)
          { int32_t ich = ch[c];
            if ((ich >= 0) && (ich < NCI) && (isnan(vmin->e[c]) || (isnan(vmin->e[c]))))
              { float_image_update_sample_range(fim, ich, &smin, &smax); }
          }
        demand(vmin < vmax, "cannot choose default scaling range");
        /* Set alll unspecified ranges: */
        for (int32_t c = 0; c < nc; c++)
          { if (isnan(vmin->e[c])) { vmin->e[c] = smin; }
            if (isnan(vmax->e[c])) { vmax->e[c] = smax; }
          }
      }
    else
      { for (int32_t c = 0; c < nc; c++)
          { int32_t ich = ch[c];
            if ((ich >= 0) && (ich < NCI) && (isnan(vmin->e[c]) || (isnan(vmin->e[c]))))
              { smin = +INF; smax = -INF;
                float_image_update_sample_range(fim, ich, &smin, &smax);
                demand(vmin < vmax, "cannot choose default scaling range");
                if (isnan(vmin->e[c])) { vmin->e[c] = smin; }
                if (isnan(vmax->e[c])) { vmax->e[c] = smax; }
              }
          }
       }
  }

float_image_t *ftp_read_fni_image(FILE *rd, int32_t *NC, int32_t *NX, int32_t *NY)
  { float_image_t *fim = float_image_read(rd);
    (*NC) = (int32_t)fim->sz[0];
    (*NX) = (int32_t)fim->sz[1];
    (*NY) = (int32_t)fim->sz[2];
    return fim;
  }

void ftp_write_png_image
  ( FILE *wr, 
    float_image_t *fim, 
    bool_t isMask, 
    int32_t nc, 
    int32_t ch[], 
    double lo[], 
    double hi[],
     uint16_t maxval, 
    float nanval, 
    bool_t verbose
  )
  { int32_t NCI = (int32_t)fim->sz[0];
    /* Replace {NAN}s: */
    int32_t c;
    for (c = 0; c < nc; c++)
      { int32_t ich = ch[c];
        if ((ich >= 0) && (ich < NCI))
          { float v = nanval;
            if (isnan(v)) { v = (float)(lo[c] + hi[c])/2; }
            float_image_replace_nan_samples(fim, ich, v);
          }
      }
    
    /* !!!  {float_image_to_uint16_image} should take a gamma parameter !!! */
    double gamma = 1.0;

    /* Convert {fim} to PNG image: */
    bool_t yup = TRUE; /* Y axis up. */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, nc, lo, hi, ch, maxval, yup, verbose);
    /* Write to disk: */
    uint16_image_write_png_file(wr, pim, gamma, verbose);
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
    if (argparser_keyword_present(pp, "-channel"))
      { o->channels = int32_vec_new(1);
        o->channels.e[0] = (int32_t)argparser_get_next_int(pp, INT32_MIN, INT32_MAX);
      }
    else 
      { o->channels = ftp_parse_int32_vec(pp, "-channels"); }

    /* Parse input range specs: */
    o->min = ftp_parse_double_vec(pp, "-min");
    o->max = ftp_parse_double_vec(pp, "-max");

    /* Uniform scaling option: */
    if (argparser_keyword_present(pp, "-uniform"))
      { o->uniform = argparser_get_next_bool(pp); }
    else
      { o->uniform = TRUE; }

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

    if (argparser_keyword_present(pp, "-verbose"))
      { o->verbose = argparser_get_next_bool(pp); }
    else
      { o->verbose = FALSE; }

    argparser_finish(pp);

    return o;
  }

double_vec_t ftp_parse_double_vec(argparser_t *pp, char *key)
  { double_vec_t vex = double_vec_new(4);
    int32_t nv = 0;
    if (argparser_keyword_present(pp, key))
      { while (argparser_next_is_non_keyword(pp))
          { vex.e[nv] = argparser_get_next_double(pp, -INF, +INF);
            nv++;
          }
        if (nv == 0) { argparser_error(pp, "should specify at least one value"); }
      }
    double_vec_trim(&(vex), nv);
    return vex;
  }

int32_vec_t ftp_parse_int32_vec(argparser_t *pp, char *key)
  { int32_vec_t vex = int32_vec_new(4);
    int32_t nv = 0;
    if (argparser_keyword_present(pp, key))
      { while (argparser_next_is_non_keyword(pp))
          { vex.e[nv] = (int32_t)argparser_get_next_int(pp, INT32_MIN, INT32_MAX);
            nv++;
          }
        if (nv == 0) { argparser_error(pp, "should specify at least one value"); }
      }
    int32_vec_trim(&(vex), nv);
    return vex;
  }
