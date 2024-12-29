#define PROG_NAME "pnmwiggle"
#define PROG_DESC "mix uniform random noise into a PBM/PGM/PPM file"
#define PROG_VERS "1.0"
/* Last edited on 2024-12-21 11:59:23 by stolfi */

#define PROG_COPYRIGHT "Copyright © 1996 by the State University of Campinas (UNICAMP)"

/* !!! Handle bad pixel values as in {pgmwfilter}/{pnmwfilter}. !!! */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [ -seed {SEED} ] {AMOUNT} \\\n" \
  "  [ < ] {INFILE}.pgm > {OUTFILE}.pgm\n" \
  "The AMOUNT should be a number in [0.0 _ 1.0].\n" \
  "Value 0.0 means no change, 1.0 means pure uniform noise in [0 .. maxval]."

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads a portable graymap as input, mixes" \
  " it with a specified amount of white noise, and writes it out.\n" \
  "\n" \
  "  The fraction of noise in the output is the {AMOUNT} parameter" \
  " (a real number between 0 and 1). More precisely, the output" \
  " image is {1-AMOUNT} times the input, plus {AMOUNT} times a random" \
  " noise uniformly distributed between 0 and the input's {maxval}.\n" \
  "\n" \
  "  The output image will have the same number of channels, size and {maxval} as" \
  " the input (except that the output {maxval} will be set to 255, if" \
  " if the input {maxval} is less than that).\n" \
  "\n" \
  "  The noise is derived from the random(3C) generator.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -seed {SEED}\n" \
  "    This option causes the random number generator to" \
  " be initialized with the specified seed. The default" \
  " is a fixed internal seed.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgmnoise(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  The program was created on 1996-11-21 by J. Stolfi, UNICAMP as \"pgmwiggle\".\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2009-09-29 J.Stolfi, IC-UNICAMP: general PNM input, renamed \"pnmwiggle\".\n" \
  "  2006-11-20 J.Stolfi, IC-UNICAMP: rewrote using Stolfi's PNM library.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <bool.h> 
#include <jsrandom.h> 
#include <jspnm.h> 
#include <jsfile.h> 
#include <uint16_image.h> 
#include <argparser.h> 

typedef struct options_t 
  { char *imgname;
    double amount;
    int seed;
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char* argv[]);
options_t *parse_options(int argc, char **argv);

/* ROUTINES */
 
int main(int argc, char* argv[])
  {
    options_t *o = parse_options(argc, argv);

    /* Image dimensions: */
    int chns, cols, rows;

    /* Open input image file: */
    uint16_t imaxval;
    bool_t iraw, ibits;
    pnm_format_t iformat;
    FILE *infile = open_read(o->imgname, TRUE);
    pnm_read_header(infile, &cols, &rows, &chns, &imaxval, &iraw, &ibits, &iformat);
    uint16_t *ipix = uint16_image_alloc_pixel_row(cols, chns);
    
    /* Choose output format and write output image header: */
    uint16_t omaxval = (imaxval < 255 ? 255 : imaxval);
    bool_t oraw, obits;
    pnm_format_t oformat;
    pnm_choose_output_format(omaxval, chns, FALSE, &oformat, &oraw, &obits);
    pnm_write_header(stdout, cols, rows, omaxval, oformat);
    
    uint16_t *opix = uint16_image_alloc_pixel_row(cols, chns);
    
    int x, y;
    double amount = o->amount, tnuoma = 1 - amount;
    srandom(o->seed);
    for (y = rows-1; y >= 0; y--)
      { pnm_read_pixels(infile, ipix, cols, chns, imaxval, iraw, ibits);
        for (x = 0; x < chns*cols; x++)
          { double iv = (ipix[x] + 0.5)/((double)imaxval + 1.0);
            double ov = tnuoma*iv + amount*drandom();
            bool_t isMask = FALSE; /* Quantize as a smooth-valued image. */
            opix[x] = pnm_quantize(ov, omaxval, isMask, PNM_NO_BADVAL);
          }
        pnm_write_pixels(stdout, opix, cols, chns, omaxval, oraw, obits);
      }
    fclose(infile);
    fclose(stdout);
    return(0);
  } 

options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    o->imgname = NULL;
    
    
    if (argparser_keyword_present(pp, "-seed"))
      { o->seed = (int)argparser_get_next_int(pp, 1, INT_MAX); }
    else
      { o->seed = 46150417; }


    argparser_skip_parsed(pp);
    o->amount = argparser_get_next_double(pp, 0.0, 1.0);
    
    if (argparser_next(pp) != NULL)
      { o->imgname = argparser_get_next(pp); }
    else
      { o->imgname = "-"; }
    
    argparser_finish(pp);
    return o;
  }
