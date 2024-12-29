#define PROG_NAME "pikewhitemask"
#define PROG_DESC "build an AVC Pike F-100 white mask"
#define PROG_VERS "1.0"

#define pikewhitemask_C_COPYRIGHT "Copyright © 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2024-12-21 12:00:33 by stolfi */

#define PROG_HELP \
  PROG_NAME  " \\\n" \
  "  [ -equalize {COL} {GAIN0} {GAIN1} ] \\\n" \
  "  [ -verbose ] \\\n" \
  "  [ < ] {RAW16FILE} > {PPMFILE}"

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
  "  pgm(5), ppm(5), pnm(5), pgmselect(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in jun/2010 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "  Info and experiments by Rafael Saracchini, IC-UNICAMP/MVL-UWE." \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Option bla bla added by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pikewhitemask_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program reads a 1000x1000 raw image file produced by the AVC Pike F-100" \
  " camera (which uses a Kodak KAI-1020 sensor), in {RAW16} mode.  The image" \
  " must be a photo of a smooth white surface taken with the" \
  " lens out of focus.  Outputs a white field mask suitable for pixel gain correction.\n" \
  "\n" \
  "  The mask is written as a PGM image with {maxval = 65535}, and" \
  "\n" \
  " samples in the range {maxval/2..maxval}.  "
  
#define PROG_INFO_OPTS \
  "  -equalize {COL} {GAIN0} {GAIN1}\n" \
  "    This option is meant to compensate for inproper gains in" \
  " the two output amplifiers of the KAI-1020 chip.  It specifies floating" \
  " point factors {GAIN0,GAIN1} to be multiplied by the sample" \
  " values in the left half of the image (columns {0..COL-1}) and on the" \
  " right half (columns {COL..999}), respectively.  (According to the" \
  " KAI-1020 documentation the column {COL} should be 500 but for" \
  " some reason it is 502 in the RAW16 files produced by the camera.)  By" \
  " default, no output amp equalization is done.\n" \
  "\n" \
  "  -verbose\n" \
  "    Produces diagnostic output."

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <frgb_ops.h>
#include <affirm.h>
#include <argparser.h>
#include <uint16_image.h>
#include <jspnm.h>

#include <uint16_image_Pike_F100.h>

/* PROTOTYPES */

typedef struct options_t
  { /* Files: */
    char *input_file_name;   /* Input filename, or "-" if not given. */
    /* Filters: */
    int equalize_col;       /* Break column for output amp equalization; or -1 if none. */
    double equalize_gain0;  /* Left-half factor for output amp equalization. */
    double equalize_gain1;  /* Right-half factor for output amp equalization. */
    /* Diagnostics: */
    bool_t verbose;
  } options_t;
  /* Command line arguments. */

int main(int argc, char* argv[]);

options_t *pikewhitemask_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

/* REAL CODE */

int main(int argc, char** argv)
  {
    /* Parse command line arguments. */
    options_t *o = pikewhitemask_parse_options(argc, argv);

    /* Get input image: */
    uint16_image_t *img = uint16_image_Pike_F100_read(o->input_file_name, o->verbose);
    assert(img->chns == 1);
    
    if (o->equalize_col >= 0)
      { demand(o->equalize_col < img->cols, "invalid column in \"-equalize\""); 
        uint16_image_Pike_F100_output_amp_balance
          ( img, o->equalize_col, o->equalize_gain0, o->equalize_gain1, o->verbose );
      }
      
    uint16_image_t *omg = pnm_Pike_F100_bayer_channel_white_mask(img, o->verbose);
      
    /* Write output image: */
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(stdout, omg, forceplain, o->verbose);
    
    /* All done: */
    if (o->verbose) { fprintf(stderr, "done.\n"); }
    exit(0);
  }

options_t *pikewhitemask_parse_options(int argc, char **argv)
  {
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    if (argparser_keyword_present(pp, "-equalize"))
      { o->equalize_col = argparser_get_next_int(pp, -1, 999);
        o->equalize_gain0 = argparser_get_next_double(pp, 0.0, 100.0);
        o->equalize_gain1 = argparser_get_next_double(pp, 0.0, 100.0);
      }
    else 
      { o->equalize_col = -1;
        o->equalize_gain0 = NAN;
        o->equalize_gain1 = NAN;
      }
    
    o->verbose = argparser_keyword_present(pp, "-verbose");
    
    /* Get optional input file name: */
    argparser_skip_parsed(pp);
    if (argparser_next(pp) != NULL) 
      { o->input_file_name = argparser_get_next(pp); }
    else
      { o->input_file_name = "-"; }

    /* Check for extraneous arguments: */
    argparser_finish(pp);
    
    return o;
  }
