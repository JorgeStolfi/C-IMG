#define PROG_NAME "pnmxarith"
#define PROG_DESC "perform arithmetic on two PBM/PGM/PPM image files"
#define PROG_VERS "2.0"

/* Copyright © 1989, 1991 by Jef Poskanzer.
** See end of file for full copyright and (no)warranty notice.
** Last edited on 2017-10-26 18:48:56 by stolfilocal
*/

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    { -add | -subtract | -difference | -diff2 | -maximum | -minimum \\\n" \
  "    | -multiply | -divide | -ratio \\\n" \
  "    | -mix {ALPHA} {BETA} \\\n" \
  "    } \\\n" \
  "    [ -scale {SCALE} ] [ -offset {OFFSET} ] \\\n" \
  "    " argparser_help_info_HELP " \\\n" \
  "    {PNMFILE1} {PNMFILE2} \\\n" \
  "    > {PNMFILE_OUT} \\\n"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "   Reads two PBM/PGM/PPM image files as input.  Performs" \
  " the specified arithmetic operation on corresponding samples," \
  " and produces a portable anymap as output.\n" \
  "\n" \
  "   The two input images must be the same width and height." \
  " Either operand gets promoted to PGM or PPM and deepened" \
  " as necessary to match the other's type and {maxval} (max sample value).\n" \
  "\n" \
  "   If either {PNMFILE1} or {PNMFILE2} is \"-\", the image" \
  " is read from standard input.   If the arguments {PNMFILE1} and {PNMFILE2}" \
  " are identical (in particular, if both are \"-\"), the image file" \
  " is read only once, and used for both operands.\n" \
  "\n" \
  "   The arithmetic is performed between corresponding pixels in the two" \
  " images, channel by channel. The input pixel values are" \
  " divided by the input {maxval}" \
  " before the operation, so that they range in [0..1]. The result is" \
  " multiplied by the given scale factor {SCALE} (default 1.0)" \
  " and added to the {OFFSET} (default 0)." \
  " The result is then scaled by the output {maxval}" \
  " rounded and clipped to the range [0..{maxval}].\n" \
  "\n" \
  "   This program is similar to \"pnmarith\", with added operations.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -scale {SCALE}\n" \
  "    Specifies an extra scale factor for the output" \
  " samples. The default is 1.0.\n" \
  "\n" \
  "  -offset {OFFSET}\n" \
  "    Specifies an offset to be added to every output" \
  " sample, after the \"-scale\" option but before scaling" \
  " by the output {maxval}. The default is 0.0.\n" \
  "\n" \
  "OPERATIONS\n" \
  "  -difference\n" \
  "    Computes the absolute value of {PNMFILE1} minus" \
  " {PNMFILE2}.\n" \
  "\n" \
  "  -diff2\n" \
  "    Computes the square of the difference {PNMFILE1} minus" \
  " {PNMFILE2}.\n" \
  "\n" \
  "  -subtract\n" \
  "    Computes the signed difference {PNMFILE1} minus" \
  " {PNMFILE2}. (The \"-offset\" parameter should then be" \
  " used to avoid loss of negative values.)\n" \
  "\n" \
  "  -ratio\n" \
  "    Computes {PNMFILE1} minus {PNMFILE2}" \
  " divided by the sum {PNMFILE1+PNMFILE2}.\n" \
  "\n" \
  "  -mix {ALPHA} {BETA}\n" \
  "    Computes the linear combination of {PNMFILE1} and {PNMFILE2} with weights" \
  " {Alpha} and {BETA}, respectively.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmarith(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in 1989-1991 by Jef Poskanzer as \"pnmarith\".\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  nov/2006: Rewritten to use lean PBM libs, argparser, etc. By Jorge Stolfi. \n" \
  "\n" \
  "  apr/2006: If files are the same, read only one. By Jorge Stolfi.  \n" \
  "\n" \
  "  nov/2001: Options \"-ratio\", \"-mix\" added by Jorge Stolfi.\n" \
  "\n" \
  "  (unknown date): Options \"-div\", \"-max\", \"-min\", \"-scale\", \"-offset\"" \
  "    added by Jorge Stolfi.\n" \
  "\n" \
  "  (unknown date): Slightly modified by Marcel Wijkstra <wijkstra@fwi.uva.nl>.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 1989, 1991 by Jef Poskanzer.\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <values.h>

#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <argparser.h>
#include <bool.h>

#define MAXCHANNELS 3

#define INF INFINITY

/* COMMAND-LINE OPTIONS */

typedef struct rgb_pixel_t { uint16_t c[3]; } rgb_pixel_t; 

typedef struct options_t
  { char function;
    double scale;
    double offset;
    double alpha;
    double beta;
    char *fname1;
    char *fname2;
  } options_t;

/* INTERNAL PROTOTYPES */

options_t *parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int main(int argc,char** argv);

/* Compute "a/b" scaled by "maxv". */
#define XELDIV(a,b,maxv) \
  ((xelval)((((double)(a))/((double)(b)))*(maxv) + 0.5))

extern double strtod(const char *, char **);

/* IMPLEMENTATIONS */

int main(int argc, char** argv)
  {
    options_t *o = parse_options(argc, argv);

    /* Parameters of image file 1: */
    FILE* ifp1; 
    int rows1, cols1, chns1;
    uint16_t maxval1;
    bool_t raw1, bits1;
    pnm_format_t format1;
    uint16_t* smp1;
    double scale1;
    
    /* Parameters of image file 2: */
    FILE* ifp2;
    int rows2, cols2, chns2;
    uint16_t maxval2;
    bool_t raw2, bits2;
    pnm_format_t format2;
    uint16_t* smp2;
    double scale2;
   
    /* Open files and allocate row buffer (only once if they are the same): */
    ifp1 = open_read(o->fname1, FALSE);
    pnm_read_header(ifp1, &cols1, &rows1, &chns1, &maxval1, &raw1, &bits1, &format1);
    smp1 = uint16_image_alloc_pixel_row(cols1, chns1);
    scale1 = (double)maxval1;

    if (strcmp(o->fname1,o->fname2) == 0)
      { ifp2 = NULL; 
        cols2 = cols1;
        rows2 = rows1;
        chns2 = chns1;
        maxval2 = maxval1;
        format2 = format1;
        smp2 = smp1;
      }
    else
      { ifp2 = open_read(o->fname2, FALSE);
        pnm_read_header(ifp2, &cols2, &rows2, &chns2, &maxval2, &raw2, &bits2, &format2);
        if (cols2 != cols1 || rows2 != rows1)
          { pnm_error("the two anymaps must be the same width and height"); }
        smp2 = uint16_image_alloc_pixel_row(cols2, chns2);
      }
    scale2 = (double)maxval2;

    /* Select output maxval and format: */
    FILE *ifp3;
    int cols3 = cols1;
    int rows3 = rows1;
    int chns3 = (chns1 > chns2 ? chns1 : chns2);
    uint16_t maxval3 = (maxval1 > maxval2 ? maxval1 : maxval2);
    bool_t raw3, bits3;
    pnm_format_t format3;
    uint16_t* smp3;
    double scale3;
    
    ifp3 = stdout;
    pnm_choose_output_format(maxval3, chns3, (! (raw1||raw2)), &format3, &raw3, &bits3);
    pnm_write_header(ifp3, cols3, rows3, maxval3, format3);
    smp3 = uint16_image_alloc_pixel_row(cols3, chns3);
    scale3 = (double)maxval3;
    
    if ((chns1 != chns3) || (chns2 != chns3))
      { assert(chns1 <= chns3);
        if (chns1 < chns3) { pnm_message("promoting first file to %d channels", chns3); }
        assert(chns2 <= chns3);
        if (chns2 < chns3) { pnm_message("promoting second file to %d channels", chns3); }
      }

    /* Promoted pixels (all indexed [0..cnhs3-1]): */
    double v1[MAXCHANNELS]; /* Pixel samples from image 1. */
    double v2[MAXCHANNELS]; /* Pixel samples from image 2. */
    double v3[MAXCHANNELS]; /* Pixel samples from result image. */
    
    /* Parameters for samples scaled to [0_1]: */
    double falpha = o->alpha;
    double fbeta = o->beta;
    double fscale = o->scale;
    double foffset = o->offset;
    
    int row, col, c;
    uint16_t *x1P;
    uint16_t *x2P;
    uint16_t *x3P;
    for (row = 0; row < rows3; ++row)
      {
        pnm_read_pixels(ifp1, smp1, cols1, chns1, maxval1, raw1, bits1);
        if (ifp2 != NULL)
          { pnm_read_pixels(ifp2, smp2, cols2, chns2, maxval2, raw2, bits2); }

        for (col = 0, x1P = smp1, x2P = smp2, x3P = smp3; col < cols1; ++col)
          {
            /* Get samples from image 1 and promote to {chns3} channels: */
            for (c = 0; c < chns1; c++) { v1[c] = (*x1P)/scale1; x1P++; }
            for (c = chns1; c < chns3; c++) { v1[c] = v1[0]; }

            /* Get samples from image 2 and promote to {chns3} channels: */
            for (c = 0; c < chns2; c++) { v2[c] = (*x2P)/scale2; x2P++; }
            for (c = chns2; c < chns3; c++) { v2[c] = v2[0]; }
            
            /* Apply operation: */
            switch (o->function)
              {
                case '+':
                  for (c = 0; c < chns3; c++) 
                    { v3[c] = v1[c] + v2[c]; }
                  break;

                case '-':
                  for (c = 0; c < chns3; c++) 
                    { v3[c] = v1[c] - v2[c]; }
                  break;

                case 'M':  /* alpha-beta mix */
                  for (c = 0; c < chns3; c++) 
                    { v3[c] = falpha*v1[c] + fbeta*v2[c]; }
                  break;

                case '*':
                  for (c = 0; c < chns3; c++) 
                    { v3[c] = v1[c] * v2[c]; }
                  break;

                case '/':
                  for (c = 0; c < chns3; c++) 
                    { v3[c] = (v2[c] == 0 ? NAN : v1[c] / v2[c]); }
                  break;

                case 'D': /* Absolute difference: */
                  for (c = 0; c < chns3; c++) 
                    { v3[c] = fabs(v1[c] - v2[c]); }
                  break;

                case '2': /* Diff squared: */
                  for (c = 0; c < chns3; c++) 
                    { double d = v1[c] - v2[c]; v3[c] = d*d; }
                  break;

                case 'R':
                  for (c = 0; c < chns3; c++) 
                    { double d = v1[c] + v2[c];
                      v3[c] = (d == 0 ? 0.5 : v1[c] / d);
                    }
                  break;

                case 'X':
                  for (c = 0; c < chns3; c++) 
                    { v3[c] = (v1[c] > v2[c] ? v1[c] : v2[c]); }
                  break;

                case 'N':
                  for (c = 0; c < chns3; c++) 
                    { v3[c] = (v1[c] < v2[c] ? v1[c] : v2[c]); }
                  break;

                default:
                  assert(FALSE);
              }

            /* Apply user scale-factor and offset: */
            for (c = 0; c < chns3; c++) { v3[c] = v3[c] * fscale + foffset; }
            
            /* Clip to [0_1] range, quantize and store: */
            for (c = 0; c < chns3; c++) 
              { if (v3[c] <= 0) 
                  { (*x3P) = 0; }
                else if (v3[c] >= 1.0) 
                  { (*x3P) = maxval3; }
                else
                  { int iv = (int)floor(v3[c] * scale3 + 0.5);
                    assert(iv <= maxval3);
                    assert(iv >= 0);
                    (*x3P) = (uint16_t)iv;
                  }
                x3P++;
              }
          }
        pnm_write_pixels(stdout, smp3, cols3, chns3, maxval3, raw3, bits3);
      }

    fflush(ifp3);
    if ((ifp1 != NULL) && (ifp1 != stdin)) { fclose(ifp1); }
    if ((ifp2 != NULL) && (ifp1 != stdin)) { fclose(ifp2); }
    if (ifp3 != NULL) { fclose(ifp3); }

    exit(0);
  }

options_t *parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    /* Parse keyword-based arguments: */
    if (argparser_keyword_present(pp, "-add"))
      { o->function = '+'; } 
    else if (argparser_keyword_present(pp, "-subtract"))
      { o->function = '-'; } 
    else if (argparser_keyword_present(pp, "-multiply"))
      { o->function = '*'; } 
    else if (argparser_keyword_present(pp, "-divide"))
      { o->function = '/'; } 
    else if (argparser_keyword_present(pp, "-difference"))
      { o->function = 'D'; } 
    else if (argparser_keyword_present(pp, "-diff2"))
      { o->function = '2'; } 
    else if (argparser_keyword_present(pp, "-ratio"))
      { o->function = 'R'; } 
    else if (argparser_keyword_present(pp, "-maximum"))
      { o->function = 'X'; } 
    else if (argparser_keyword_present(pp, "-minimum"))
      { o->function = 'N'; } 
    else if (argparser_keyword_present(pp, "-mix"))
      { o->function = 'M'; 
        o->alpha = argparser_get_next_double(pp, -MAXDOUBLE, +MAXDOUBLE);
        o->beta  = argparser_get_next_double(pp, -MAXDOUBLE, +MAXDOUBLE);
      }
    else
      { argparser_error(pp, "no operation was specified"); }
    
    if (argparser_keyword_present(pp, "-scale"))
      { o->scale = argparser_get_next_double(pp, -MAXDOUBLE, +MAXDOUBLE); }
    else
      { o->scale = 1.0; }
      
    if (argparser_keyword_present(pp, "-offset"))
      { o->offset = argparser_get_next_double(pp, -MAXDOUBLE, +MAXDOUBLE); }
    else
      { o->offset = 0.0; }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* First operand: */
    if (argparser_next(pp) != NULL)
      { o->fname1 = argparser_get_next(pp); }
    else
      { argparser_error(pp, "missing first operand"); }

    if (argparser_next(pp) != NULL)
      { o->fname2 = argparser_get_next(pp); }
    else
      { argparser_error(pp, "missing second operand"); }

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }
