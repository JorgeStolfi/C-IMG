#define PROG_NAME "pnm_to_fni"
#define PROG_DESC "convert a PGM or PPM image file to float-valued FNI file"
#define PROG_VERS "1.0"

#define PROG_C_COPYRIGHT "Copyright © 2005  State University of Campinas (UNICAMP)."

/* Last edited on 2017-06-22 19:05:50 by stolfilocal */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ " pst_scaling_parse_min_any_HELP " ] \\\n" \
  "    [ " pst_scaling_parse_max_any_HELP " ] \\\n" \
  "    [ " pst_scaling_parse_center_any_HELP " ] \\\n" \
  "    [ " pst_scaling_parse_width_any_HELP " ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    [ -isMask {ISMASK} ] \\\n" \
  "    < {PNM_FILE} > {FNI_FILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads from standard input a Portable Bitmap file, either color (PPM)" \
  " or grayscale (PGM), and writes it out as a float-valued image" \
  " (FNI) file with one or three channels, respectively.  See " \
  " {float_image.h} for a description of the FNI format.\n" \
  "\n" \
  "  Let {MAXVAL} be the maximum input sample value, as specified" \
  " in the input file's header. The input samples are first converted" \
  " from the range {0..MAXVAL} to the float interval {[0_1]}, as determined" \
  " by the \"-isMask\" option.  Then the range {[0_1]} is mapped" \
  " affinely (linearly) to some real range" \
  " [{VMIN} _ {VMAX}], specified by the other options.  The" \
  " encoding is essentially linear, without any gamma mapping.\n" \
  "\n" \
  "SPECIFYING THE SAMPLE SCALING\n" \
  "  The output scaling range [{VMIN} _ {VMAX}] is specified by" \
  " the options " pst_scaling_option_list_INFO ".  Exactly two" \
  " of the four values {VMIN}, {VMAX}, {VCTR} and {VWID} must" \
  " be specified.  " pst_scaling_complete_params_INFO "\n" \
  "\n" \
  "OPTIONS\n" \
  "  " pst_scaling_parse_min_one_HELP "\n" \
  "  " pst_scaling_parse_min_RGB_HELP "\n" \
  "    " pst_scaling_parse_min_INFO \
  " of the output scaling range.  " \
  pst_double_vec_spec_den_INFO "  " \
  pst_scaling_num_values_INFO "\n" \
  "\n" \
  "  " pst_scaling_parse_max_one_HELP "\n" \
  "  " pst_scaling_parse_max_RGB_HELP "\n" \
  "    " pst_scaling_parse_max_INFO \
  " of the output scaling range.  " \
  pst_double_vec_spec_den_INFO "  " \
  pst_scaling_num_values_INFO "\n" \
  "\n" \
  "  " pst_scaling_parse_center_one_HELP "\n" \
  "  " pst_scaling_parse_center_RGB_HELP "\n" \
  "    " pst_scaling_parse_center_INFO \
  " of the output scaling range.  " \
  pst_double_vec_spec_den_INFO "  " \
  pst_scaling_num_values_INFO "\n" \
  "\n" \
  "  " pst_scaling_parse_width_one_HELP "\n" \
  "  " pst_scaling_parse_width_RGB_HELP "\n" \
  "    " pst_scaling_parse_width_INFO \
  " of the output scaling range.  " \
  pst_double_vec_spec_den_INFO "  " \
  pst_scaling_num_values_INFO "\n" \
  "\n" \
  "  -isMask {ISMASK}\n" \
  "    This optional Boolean argument specifies the interpretation" \
  " of integer sample values in input file, specifically" \
  " how the integer values {0..MAXVAL} are mapped to the range {[0_1]} before" \
  " being affinely mapped to {[VMIN_VMAX]}.  If {ISMASK} is true (\"T\" or 1)," \
  " " sample_conv_0_1_isMask_true_INFO "  If {ISMASK} is false (\"F\" or 0)," \
  " " sample_conv_0_1_isMask_false_INFO "  The default is \"-isMask F\".\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  fni_to_pnm(1), float_image.h.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2005-12-10 by Jorge Stolfi, UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2005-12-10 created.\n" \
  "  2010-08-14 added the \"-isMask\" flag.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <assert.h>

#include <argparser.h>
#include <float_image.h>
#include <float_image_from_uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <sample_conv.h>
#include <bool.h>

#include <pst_scaling.h>
#include <pst_basic.h>

#define INF INFINITY

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { int NC;              /* Number of channels expected in input, or -1 if unknown. */
    bool_t isMask;       /* Interpretation of integer sample values. */
    /* Output channel data: */
    double_vec_t min;    /* Low endpoint of scaling range; empty vec if not given. */
    double_vec_t max;    /* High endpoint of scaling range; empty vec if not given. */
    double_vec_t ctr;    /* Center of scaling range; empty vec if not given. */
    double_vec_t wid;    /* Width of scaling range; empty vec if not given. */
  } options_t;

/* PROTOTYPES */

void write_uint16_image_as_float_image(uint16_image_t *pim, bool_t isMask, double min, double max);
  /* Writes the grayscale PNM image {pim} to stdout as a
    single-channel float image, mapping pixel values from {[0 ..
    maxval]} to {[min _ max]}. */ 

void write_ppm_image_as_float_image(uint16_image_t *pim, bool_t isMask, double_vec_t *min, double_vec_t *max);
  /* Writes the color PNM image {pim} to stdout as a three-channel
    float image, mapping pixel values in channel {c} from {[0 ..
    maxval]} to {[min[c] _ max[c]]}. */ 

int main(int argc, char** argv);

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

/* IMPLEMENTATIONS */

#define ptf_debug TRUE

int main(int argc, char** argv)
  { 
    /* Parse command line arguments: */
    options_t *o = parse_options(argc, argv);
    
    /* Read input image: */
    uint16_image_t *pim = uint16_image_read_pnm_file(stdin);

    /* Check the number of channels {NC}: */
    int NC = pim->chns;
    demand((NC == 1) || (NC == 3), "input image must be PGM or PPM"); 
    demand
      ( (o->NC == -1) || (o->NC == 1) || (o->NC == NC), 
        "image channels are inconsistent with args"
      ); 
    
    /* Complete the output scaling parameters: */
    pst_scaling_fix_params
      ( NC, FALSE,
        &(o->min), &(o->max), &(o->ctr), &(o->wid),
        NULL, NULL
      );

    /* Convert the PGM or PPM image {pim} to FNI and write the result: */
    if (NC == 1)
      { write_uint16_image_as_float_image(pim, o->isMask, o->min.e[0], o->max.e[0]); }
    else
      { write_ppm_image_as_float_image(pim, o->isMask, &(o->min), &(o->max)); }
 
    return 0;
  }

void write_uint16_image_as_float_image(uint16_image_t *pim, bool_t isMask, double min, double max)
  { demand(pim->chns == 1, "input image must be PGM"); 
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, &min, &max, TRUE, TRUE);
    float_image_write(stdout, fim);
    float_image_free(fim);
  }

void write_ppm_image_as_float_image(uint16_image_t *pim, bool_t isMask, double_vec_t *min, double_vec_t *max)
  { demand(pim->chns == 3, "input image must be PPM"); 
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, min->e, max->e, TRUE, TRUE);
    float_image_write(stdout, fim);
    float_image_free(fim);
  }

options_t *parse_options(int argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    o->NC = -1; /* Number of channels deduced from the command line args. */
      
    /* Parse the integer decoding options: */
    if (argparser_keyword_present(pp, "-isMask"))
      { o->isMask = argparser_get_next_bool(pp); }
    else
      { o->isMask = FALSE; }
    
    /* Parse the output range specs: */
    o->min = pst_scaling_parse_range_option(pp, "-min",    &(o->NC));
    o->max = pst_scaling_parse_range_option(pp, "-max",    &(o->NC));
    o->ctr = pst_scaling_parse_range_option(pp, "-center", &(o->NC));
    o->wid = pst_scaling_parse_range_option(pp, "-width",  &(o->NC));

    argparser_finish(pp);
    
    return o;
  }
