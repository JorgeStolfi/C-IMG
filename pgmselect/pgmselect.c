#define PROG_NAME "pgmselect"
#define PROG_DESC "make a mask for pixels in a given value range"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-21 12:00:44 by stolfi */

/* Copyright © 1996 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -range {LO} {HI} \\\n" \
  "    < INFILE.pgm \\\n" \
  "    > OUTFILE.pgm"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads a portable graymap (PGM) image file as input, and" \
  " writes to stdout a greymap of the same size and {maxval} which" \
  " is equal to {maxval} at every pixel where the" \
  " corresponding input image pixel lies in a specified" \
  " range, and 0 elsewhere.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -range {LO} {HI}\n" \
  "    This mandatory argument specifies that the" \
  " range of input sample values to select" \
  " is {LO..HI}, including both.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgm(5), and /Dr. Strangelove/ by Stanley Kubrick(1), if you haven't seen it already.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 1996-11-21 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-03-03 J.Stolfi: Embedded the manpage (this {PROG_INFO} string).\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 1996 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include <bool.h> 
#include <jspnm.h> 
#include <jsfile.h> 
#include <uint16_image.h> 
#include <argparser.h> 

typedef struct options_t 
  { int32_t range_LO;
    int32_t range_HI;
  } options_t;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);
options_t *get_options(int32_t argc, char **argv);

/* ROUTINES */
 
int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    /* Read input image header: */
    int32_t chns, cols, rows;
    uint16_t maxval;
    bool_t raw, bits;
    pnm_format_t format;
    pnm_read_header(stdin, &cols, &rows, &chns, &maxval, &raw, &bits, &format);
    if (chns != 1) { pnm_error("input image must be monochromatic"); }
    
    /* Write output image header: */
    pnm_write_header(stdout, cols, rows, maxval, format);
    
    uint16_t *ipix = uint16_image_alloc_pixel_row(cols, chns);
    uint16_t *opix = uint16_image_alloc_pixel_row(cols, chns);

    int32_t x, y;
    for (y = 0; y < rows; y++)
      { pnm_read_pixels(stdin, ipix, cols, chns, maxval, raw, bits);
        for (x = 0; x < cols; x++)
          { opix[x] = ((ipix[x] >= o->range_LO) && (ipix[x] <= o->range_HI) ? maxval : 0); }
        pnm_write_pixels(stdout, opix, cols, chns, maxval, raw, bits);
      } 
    
    fflush(stdout);
    return(0);
  } 

options_t *get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 

    argparser_get_keyword(pp, "-range");
    o->range_LO = (int32_t)argparser_get_next_int(pp, 1, PNM_FILE_MAX_MAXVAL);
    o->range_HI = (int32_t)argparser_get_next_int(pp, 1, PNM_FILE_MAX_MAXVAL);

    argparser_skip_parsed(pp);
    
    argparser_finish(pp);

    return o;
  }
