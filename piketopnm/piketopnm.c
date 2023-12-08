#define PROG_NAME "piketopnm"
#define PROG_DESC "convert AVC Pike F-100 RAW16 image files to PPM/PGM"
#define PROG_VERS "1.0"

#define piketopnm_C_COPYRIGHT "Copyright © 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2017-06-30 01:04:42 by stolfilocal */

#define PROG_HELP \
  PROG_NAME  " \\\n" \
  "  [ -crop {COL0} {ROW0} {NCOLS} {NROWS} ] \\\n" \
  "  [ -debayer [ pad | interpolate | squeeze ] | \\\n" \
  "    -extract { R | G0 | G1 | B } \\\n" \
  "  ] \\\n" \
  "  [ -balance " frgb_parse_color_HELP " ] \\\n" \
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
  "  2010-06-22 J. Stolfi: added \"-crop\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " piketopnm_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program reads a 1000x1000 raw image file produced by the AVC Pike F-100" \
  " camera (which uses a Kodak KAI-1020 sensor), in {RAW16} mode.  It outputs" \
  " a portable pixmap {PPMFILE}.  It optionally applies a de-Bayering filter" \
  " and other corrections.\n" \
  "\n" \
  "  If neither \"-debayer\" nor \"-extract\" arguments are given, the output image will be" \
  " a PGM image with the R, G, and B samples interleaved in the camera's 2x2" \
  " Bayer patern ((G0,R0),(B0,G1)).  "
  
#define PROG_INFO_OPTS \
  "  -equalize {COL} {GAIN0} {GAIN1}\n" \
  "    This option is meant to compensate for inproper gains in" \
  " the two output amplifiers of the KAI-1020 chip.  It specifies floating" \
  " point factors {GAIN0,GAIN1} to be multiplied by the sample" \
  " values in the left half of the image (columns {0..COL-1}) and on the" \
  " right half (columns {COL..999}), respectively.  The" \
  " column {COL} refers to the whole image, before any cropping.  (According to the" \
  " KAI-1020 documentation the column {COL} should be 500 but for" \
  " some reason it is 502 in the RAW16 files produced by the camera.)  Beware" \
  " that gains above 1 may cause some samples to overflow the maximum" \
  " output value (65535); any such samples will be set to 65535.  This option is effective" \
  " even in compbination with \"-debayer\" or \"-extract\".  By" \
  " default, no output amp equalization is done.\n" \
  "\n" \
  "  -crop {ICOL} {IROW} {NCOLS} {NROWS}\n" \
  "    If this option is present, the input image is cropped to" \
  " the rectange that starts in column {ICOL} and row {IROW}, and" \
  " has {NCOLS} columns and {NROWS} rows.  All four arguments" \
  " must be even.  The default is no cropping.\n" \
  "\n" \
  "  -debayer [ pad | interpolate | squeeze ]\n" \
  "    This optional argument requests and specifies the unraveling of the" \
  " Bayer pattern.  If present, the output will be a three channel PPM" \
  " image, and each sample of the input file will be moved to the appropriate" \
  " channel of the output.  This option is mutually exclusive with \"-extract\".\n" \
  "\n" \
  "    With the \"pad\" option, the output image will have the same size" \
  " as the input (after cropping), and the samples will retain their" \
  " original column and row indices.  Thus, each pixel will have only" \
  " one significant sample; and each block of 2 by 2 pixels" \
  " in the output will contain only one significant sample in the" \
  " red and blue channels, and two in the green channel.  The two non-significant samples of" \
  " each output pixel will be set to zero.\n" \
  "\n" \
  "    The \"interpolate\" option is simiar to \"pad\", except that the" \
  " non-significant samples are filled by bilinear interpolation of their" \
  " two or four nearest significant neighbors, instead of being set to zero.\n" \
  "\n" \
  "    With the \"squeeze\" option, the output image will have only half as many rows" \
  " and columns as the (cropped) input, and all its samples will be" \
  " significant.  Each sample in column {col} and row {row} of the cropped" \
  " input image will be copied to the appropriate channel in column" \
  " {col/2} and row{row/2} of the output, both rounded down.  The two" \
  " green Bayer channels G0 and G1 will be averaged together; the average" \
  " will be rounded down or up in checkerboard fashion.\n" \
  "\n" \
  "  -extract { R0 | G0 | G1 | B0 }\n" \
  "    If this argument is given, the output image will contain only the" \
  " specified Bayer channel extracted from the raw image; namely, only" \
  " one sample out of each instance of the 2x2 Bayer template ((G0,R0),(B0,G1)).  In particular," \
  " the G0 option select the samples with two even indices, while the G1" \
  " option selects the samples with two odd indices.  In any case the output" \
  " will be a PGM (grayscale) image with half as many columns and" \
  " rows as the cropped input image.   The red and" \
  " blue channels can be specified also as \"R\" and \"B\". This option is" \
  " mutually exclusive with \"-debayer\".\n" \
  "\n" \
  "  -balance " frgb_parse_color_HELP "\n" \
  "    This option specifies factors to be multiplied by the sample" \
  " values in the three RGB channels.  The factors are specified" \
  " as " frgb_parse_color_INFO "  This option is effective" \
  " even in compbination with \"-debayer\" or \"-extract\".  The default" \
  " is \"-balance 1 1 1\" (no color adjustment).\n" \
  "\n" \
  "  -verbose\n" \
  "    Produces diagnostic output."

#define _GNU_SOURCE
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

typedef enum
  { debayer_op_PAD,
    debayer_op_INTERPOLATE,
    debayer_op_SQUEEZE,
    debayer_op_NONE
  } debayer_op_t;
  /* De-Bayering method to use. */

typedef struct options_t
  { /* Files: */
    char *input_file_name;   /* Input filename, or "-" if not given. */
    /* Cropping: */
    int crop_icol;           /* Initial column for cropping or 0 if whole. */
    int crop_irow;           /* Initial row for cropping or 0 if whole. */
    int crop_cols;           /* Column count for cropping or -1 if whole. */
    int crop_rows;           /* Row count for cropping or -1 if whole. */
    /* Filters: */
    debayer_op_t debayer_op; /* De-Bayering method, or {debayer_op_NONE}. */
    bool_t extract;          /* If true, output only one Bayer channel. */
    int extract_col;         /* Selected Bayer template column (0 or 1). */
    int extract_row;         /* Selected Bayer template row (0 or 1). */
    frgb_t balance;          /* Color balance factors. */
    int equalize_col;        /* Break column for output amp equalization; or -1 if none. */
    double equalize_gain0;   /* Left-half factor for output amp equalization. */
    double equalize_gain1;   /* Right-half factor for output amp equalization. */
    /* Diagnostics: */
    bool_t verbose;
  } options_t;
  /* Command line arguments. */

int main(int argc, char* argv[]);

options_t *piketopnm_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

/* REAL CODE */

int main(int argc, char** argv)
  {
    /* Parse command line arguments. */
    options_t *o = piketopnm_parse_options(argc, argv);

    /* Get input image: */
    uint16_image_t *img = uint16_image_Pike_F100_read(o->input_file_name, o->verbose);
    if (o->verbose)
      { fprintf(stderr, "chns = %d cols = %5d rows = %d\n", img->chns, img->cols, img->rows); }
    assert(img->chns == 1);
    
    if (o->crop_cols > 0)
      { /* The crop rectangle must be even so as not to mess up the Bayer pattern: */
        assert((o->crop_icol % 2) == 0);
        assert((o->crop_irow % 2) == 0);
        assert((o->crop_cols % 2) == 0);
        assert((o->crop_rows % 2) == 0);
        /* Replace {img} by the cropped image: */
        uint16_image_t *omg = uint16_image_crop
          ( img, 0, img->chns, o->crop_icol, o->crop_cols, o->crop_irow, o->crop_rows );
        uint16_image_free(img);
        img = omg;
      }
    
    if (o->equalize_col >= 0)
      { int col = o->equalize_col - o->crop_icol;
        uint16_image_Pike_F100_output_amp_balance
          ( img, col, o->equalize_gain0, o->equalize_gain1, o->verbose );
      }
    
    if (o->debayer_op != debayer_op_NONE)
      { assert(! o->extract);
        /* Separate pixels into proper channels: */
        bool_t squeeze = (o->debayer_op == debayer_op_SQUEEZE);
        uint16_image_t *omg = uint16_image_Pike_F100_debayer(img, squeeze, o->verbose);
        uint16_image_free(img);
        img = omg;
        if (o->debayer_op == debayer_op_INTERPOLATE)
          { /* Interpolate missing samples: */
            uint16_image_Pike_F100_interpolate_bayer(img, o->verbose);
          }
       }

    if (o->extract)
      { assert(o->debayer_op == debayer_op_NONE);
        /* Extract the selected Bayer channel: */
        uint16_image_t *omg = uint16_image_Pike_F100_extract_bayer_channel
          ( img, o->extract_col, o->extract_row, o->verbose );
        uint16_image_free(img);
        img = omg;
      }
   
    frgb_t ones = frgb_Ones;
    if (! frgb_eq(&(o->balance), &ones))
      { /* Apply color balance factors: */
        uint16_image_Pike_F100_color_balance(img, o->balance.c, o->verbose);
      }
      
    /* Write output image: */
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(stdout, img, forceplain, o->verbose);
    
    /* All done: */
    if (o->verbose) { fprintf(stderr, "done.\n"); }
    exit(0);
  }

options_t *piketopnm_parse_options(int argc, char **argv)
  {
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    if (argparser_keyword_present(pp, "-debayer"))
      { if (argparser_keyword_present_next(pp, "pad"))
          { o->debayer_op = debayer_op_PAD; }
        else if (argparser_keyword_present_next(pp, "interpolate"))
          { o->debayer_op = debayer_op_INTERPOLATE; }
        else if (argparser_keyword_present_next(pp, "squeeze"))
          { o->debayer_op = debayer_op_SQUEEZE; }
        else
          { argparser_error(pp, "unrecognized de-Bayering method"); }
      }
    else
      { o->debayer_op = debayer_op_NONE; }
    
    if (argparser_keyword_present(pp, "-extract"))
      { if (o->debayer_op != debayer_op_NONE)
          { argparser_error(pp, "cnnot use \"-extract\" with \"-debayer\""); }
        o->extract = TRUE;
        if (argparser_keyword_present_next(pp, "G0"))
          { o->extract_col = 0; o->extract_row = 0; }
        else if (argparser_keyword_present_next(pp, "R") || argparser_keyword_present_next(pp, "R0"))
          { o->extract_col = 1; o->extract_row = 0; }
        else if (argparser_keyword_present_next(pp, "B") || argparser_keyword_present_next(pp, "B0"))
          { o->extract_col = 0; o->extract_row = 1; }
        else if (argparser_keyword_present_next(pp, "G1"))
          { o->extract_col = 1; o->extract_row = 1; }
        else 
          { argparser_error(pp, "unrecognized channel"); }
      }
    else
      { o->extract = FALSE; }
    
    if (argparser_keyword_present(pp, "-balance"))
      { o->balance = frgb_parse_color(pp); }
    else 
      { o->balance = frgb_Ones; }
    
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
    
    if (argparser_keyword_present(pp, "-crop"))
      { o->crop_icol = argparser_get_next_int(pp, 0, 999);
        o->crop_irow = argparser_get_next_int(pp, 0, 999);
        o->crop_cols = argparser_get_next_int(pp, 2, 1000 - o->crop_icol);
        o->crop_rows = argparser_get_next_int(pp, 2, 1000 - o->crop_irow);
        int parity = 
          (o->crop_icol % 2) |
          (o->crop_irow % 2) |
          (o->crop_cols % 2) |
          (o->crop_rows % 2);
        if (parity != 0) { argparser_error(pp, "crop rectangle must be even"); }
      }
    else 
      { o->crop_icol = 0;
        o->crop_irow = 0;
        o->crop_cols = 0;
        o->crop_rows = 0;
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
