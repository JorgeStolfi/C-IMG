#define PROG_NAME "pnmbend"
#define PROG_DESC "apply a general deformation to a pixmap"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-25 09:32:58 by stolfi */

/* Copyright © 2000 by the State University of Campinas (UNICAMP).
** See the copyright, authorship, and warranty notice at end of file. 
*/

/* !!! TO DO: Finish! */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [ -map {PNMFILE_DX} {PNMFILE_DY} [ -range {DMIN} {DMAX} ] ] \\\n" \
  "  [ -shift [ channel {CHN} ] {DX0} {DY0} ].. \\\n" \
  "  [<] {PNMFILE_IN} \\\n" \
  "  > {PNMFILE_OUT}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program bla bla bla" \
  " bla bla bla bla {X+Y} bla bla" \
  " bla {INFILE} bla \"foobar.ppm\" bla bla bla\n" \
  "\n" \
  "  Beware that bla bla bla BLEBBLE BLOB bla" \
  " bla bla bla bla.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -blabber {AMOUNT}\n" \
  "    Blabbers for that {AMOUNT}. May also bla" \
  " bla bla bla bla bla bla bla bla bla bla bla bla" \
  " bla bla bla bla bla bla bla.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnm(5).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in jul/2000 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Option bla bla added by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2000 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

// .TH pnmbend 1 "11 mar 2000"
// .IX pnmbend
// .SH NAME
// pnmbend - apply a general deformation to an anymap
// .SH SYNOPSIS
// .B pnmbend
// .B -map
// .I pnmfileDX pnmfileDY
// .RB [ -offset
// .IR float ]
// .RB [ -scale
// .IR float ]
// .RI [ pnmfile ]
// .PP
// .B pnmbend
// .R -shift
// .I floatDX0 floatDY0
// .RI [ 
// .I floatDX1 floatDY1 floatDX2 floatDY2
// .RI ]
// .RI [ pnmfile ]
// .SH DESCRIPTION
// Reads a portable anymap 
// .I A
// from 
// .IR pnmfile,
// or from standard input if 
// .I pnmfile
// is not specified. Outputs a new image 
// .I B
// whose pixel
// .IR B[x,y]
// is a copy if pixel 
// .IR A[x+DX[
// .IR x,y],
// .IR y+DY[
// .IR x,y] 
// ], where the displacement functions 
// .IR DX , DY 
// are specified by the other parameters.
// .PP
// In the 
// .B -map 
// variant, the displacement for each pixel is specified 
// by two portable anymaps files
// .IR pnmfileDX , pnmfileDY .
// All three images must have the same width and height;
// the map files will be promoted to pgm or ppm
// format, as necessary, so as to match the input file.
// In the case of a ppm file, the displacements for each channel of the 
// input image
// are taken from the corresponding channels of the map images.
// .PP
// The output image will have the same size, type and depth as
// the input.  Note that there may be lost or undefined pixels
// around the edges.
// .PP
// The
// .B -offset
// parameter (a pixel value, default 0) will be subtracted from the displacement map 
// pixels; and the result will be divided by the 
// .B -scale 
// parameter (another pixel value, default 1).
// .PP
// In the 
// .B -shift
// variant, the displacement is the same for all pixels.
// Following the keyword there must be 
// either two numbers, 
// which are taken to be the X and
// Y displacements for all channels; or three pairs
// of numbers, each pair being the displacement
// for the corresponding channel.
// .PP
// Note that, in either case, the displacements may be fractional.
// Non-integer pixel indices imply interpolation of neighboring pixels.
// The interpolation assumes that the image was Gaussian-filtered before
// sampling, and the displacements change very little from place to
// place.
// .PP
// All flags can be abbreviated to their shortest unique prefix.
// .SH "SEE ALSO"
// pnmrotate(1), pnmcut(1), pnmshear(1), pnmscale(1), pnmbend(1), pnm(5)
// .SH AUTHOR
// Copyright (C) 2000 by Jorge Stolfi <stolfi@dcc.unicamp.br>.

#include <stdio.h>
#include <math.h>

#include <bool.h>
#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <affirm.h>
#include <argparser.h>
#include <sample.h>

typedef struct options_t 
  { /* Command line arguments: */
    MapType maptype;
    char *imIname;  /* Input file name, or "-" for {stdin}. */
    /* If map given: */
    char *mapDXname;  /* X displacement map (NULL if not specified). */
    char *mapDYname;  /* Y displacement map (NULL if not specified). */
    double dmin;      /* Displacement for map pixels equal to 0. */
    double dmax;      /* Displacement for map pixels equal to  {maxval}. */
    /* Output image dimensions: */
    int cols;         /* X size of output image, or -1 if not specified. */
    int rows;         /* Y size of output image, or -1 if not specified. */
    /* Pixel displacement: */
    double shiftX[MAX_CHANNELS];
    double shiftY[MAX_CHANNELS];
  } options_t;

options_t *parse_options(int argc, char **argv);
void parse_shift_args(argparser_t *pp, int chns, double shiftX[], double shiftY[]);
void parse_map_args(argparser_t *pp, char **mapDXnameP, char **mapDYnameP, double *offsetP, double *scaleP);

int main(int argc, char **argv)
  {
    /* Parse command line arguments: */
    options_t *o = parse_options(argc, argv);
    
    /* Input image: */
    float_image_t *imI = read_image(o->imIname, INT_MAX);
    assert(imI != NULL);
    int chns = imI.chns;
    int icols = imI.cols;
    int irows = imI.rows;

    /* Map images: */
    float_image_t *imDX = read_image(o->mapDXname, 1);
    float_image_t *imDY = read_image(o->mapDYname, 1);
    
    /* Output image: */
    int ocols = (o->cols >= 0 ? o->cols : icols);
    int orows = (o->cols >= 0 ? o->cols : irows);
    float_image_t *imO = create_image(colsNew, rowsNew, chns);

    /* Map images: */
    if ((imDX != NULL) && (imDY != NULL))
      { demand(imDY.cols == imDY.cols, "mismatched displacement map cols");
        demand(imDY.rows == imDX.rows, "mismatched displacement map rows");
      }

    float_image_fill(imO, 0.0);

    /* Distribute input pixels oer output image: */
    int irow, icol, c;
    for (irow = 0; irow < irows; ++irow)
      { /* fprintf(stderr, "distributing row %d...\n", row); */
        for (icol = 0; icol < icols; ++icol)
          { for (c = 0; c < chns; c++)
              { double iv = (double)float_image_get_sample(imI, c, irow, icol);
                
                double dx = (double)float_image_get_sample(imDX, c, irow, icol);
                double dy = (double)float_image_get_sample(imDY, c, irow, icol);
                splat(iv, 
            if (maptype == MT_IMAGE)
              { /* Obtain the displacements */
                switch (PNM_FORMAT_TYPE(imDX.format))
                  { case PPM_TYPE: 
                      { shiftX[0] = (*pDX++); 
                        shiftX[1] = (*pDX++); 
                        shiftX[2] = (*pDX++); 
                        break;
                      }
                    default:
                      { shiftX[0] = (*pDX++);
                        shiftX[1] = shiftX[0]; 
                        shiftX[2] = shiftX[0]; 
                        break;
                      }
                   }
                switch (PNM_FORMAT_TYPE(imDY.format))
                  { case PPM_TYPE: 
                      { shiftY[0] = (*pDY++); 
                        shiftY[1] = (*pDY++); 
                        shiftY[2] = (*pDY++); 
                        break;
                      }
                    default:
                      { shiftY[0] = (*pDY++);
                        shiftY[1] = shiftY[0]; 
                        shiftY[2] = shiftY[0]; 
                        break;
                      }
                   }
               }

            /* Compute the output pixel: */
            for (r = 0; r < imO.chns; r++)
              { double x, y;
                x = (double)col + shiftX[r];
                y = (double)row + shiftY[r];
                /* fprintf(stderr, " row = %d col = %d chan = %d pO = %x", row, col, r, (unsigned int)pO); */
                sample_pixel(&imI, x, y, r, pO);
                /* fprintf(stderr, " v = %.3f", *pO); */
                pO++;
              }
          }
        float_image_write_pnm_buffer_dump_first_row(stdout, &imO, row);
        /* fprintf(stderr, "done.\n"); */
      }

    if(maptype == MT_IMAGE) { if (fpDX != std?) { fclose(fpDX); } else { fflush(fpDX); }; pm_close(fpDY); }
    if (fpI != std?) { fclose(fpI); } else { fflush(fpI); };
    if (stdout != std?) { fclose(stdout); } else { fflush(stdout); };

    return(0);
  }

void parse_shift_args(argparser_t *pp, int chns, double shiftX[], double shiftY[])
  { 
    /* Parses the arguments for the "-shift" or "-colorShift" mode. */
    /* Sets {*shiftXP}, {*shiftYP}. */
    int r;
    /* Parse the X and Y shift arguments: */
    for (r = 0; r < chns; r++)
      { shiftX[r] = argparser_get_next_double(pp, -1.0e+6, +1.0e+6);
        shiftY[r] = argparser_get_next_double(pp, -1.0e+6, +1.0e+6);
      }
    /* If "-shift", replicate the X and Y shift for all three channels: */
    for (r = chns; r < 3; r++)
      { shiftX[r] = shiftX[0];
        shiftY[r] = shiftY[0];
      }
  }

void parse_map_args
  ( argparser_t *pp, 
    char **mapDXnameP, 
    char **mapDYnameP, 
    double *offsetP, 
    double *scaleP
  )
  /* Parses the arguments for the "-mam" mode. */
  { /* Sets {*mapXnameP}, {*mapYnameP}, {*scaleP}, {*offsetP}. */

    (*mapDXnameP) = argparser_get_next(pp); 
    (*mapDYnameP) = argparser_get_next(pp); 

    (*offsetP) = 0.0;
    (*scaleP) = 1.0;
    while (TRUE)
      { if (argparser_keyword_present("-scale"))
          { (*scaleP) = argparser_get_next_double(pp, 1.0e-6, 1.0e+6); }
        else if (argparser_keyword_present("-offset"))
          { (*offsetP) = argparser_get_next_double(pp, -1.0e+6, +1.0e+6); }
        else
          { break; }
      }
  }
  
options_t *get_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");

    if (argparser_keyword_present("-shift"))
      { o->maptype = MT_FIXED;
        o->shiftchans = 1;
        /* fprintf (stderr, "single-chanel shift\n"); */
        parse_shift_args(pp, o->shiftchans, &(o->shiftX[0]), &(o->shiftY[0]));
      }
    else if (argparser_keyword_present("-colorShift"))
      { o->maptype = MT_FIXED;
        o->shiftchans = 3;
        /* fprintf (stderr, "three-chanel shift\n"); */
        parse_shift_args(pp, o->shiftchans, &(o->shiftX[0]), &(o->shiftY[0]));
      }
    else if (argparser_keyword_present("-map"))
      { o->maptype = MT_IMAGE;
        parse_map_args(pp, &(o->mapDXname), &(o->mapDYname), &(o->offset), &(o->scale));
      }
    else
      { argparser_error(pp, "you must specify \"-map\" or \"-shift\""); }

    argparser_skip_parsed(pp);

    if (argparser_next(pp) == NULL) 
      { o->file_name = "-"; }
    else
      { o->file_name = argparser_get_next(pp); }

    argparser_finish(pp);

    return o;
  }
