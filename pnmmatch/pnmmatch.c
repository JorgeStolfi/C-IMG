#define PROG_NAME "pnmmatch"
#define PROG_DESC "deform a pixmap to match another pixmap"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-25 09:28:41 by stolfi */

/* Copyright © 2000 by the State University of Campinas (UNICAMP).
** See the copyright, authorship, and warranty notice at end of file.
*/

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -step NN -maxDisp NN \\\n" \
  "  -outName XXX \\\n" \
  "  [-scale NN] [ -offset NN ] \\\n" \
  "  pnmfileA pnmfileB"
  
// .TH pgmmatch 1 "11 mar 2000"
// .IX pgmmatch
// .SH NAME
// pgmmatch - compute the displacement map that best matches two images
// .SH SYNOPSIS
// .B pnmmatch
// .RB -bumpspacing
// .IR float
// .RB -maxdisp
// .IR float
// .RB [ -scale
// .IR float ]
// .RB [ -offset
// .IR float ]
// .RB [ -depth
// .IR int ]
// .RB [ -name 
// .IR string ] 
// .I pnmfile1 pnmfile2
// .SH DESCRIPTION
// Reads two portable anymaps  
// .IR A, B
// from
// .IR pnmfile1, pnmfile2.
// Finds displacement functions
// .IR DX[x,y]
// and
// .IR DY[x,y]
// such that 
// .IR A[x+DX[x,y]
// ,
// .IR y+DY[x,y]
// ] approximates 
// .IR B[x,y]
// for all pixels 
// .IR x,y.
// Note that the displacements are geenrally fractional,
// and imply interpolation of neighboring pixels.
// .PP
// Each of the maps 
// .IR DX,DY
// is the sum of Gaussian bumps centered 
// at the corners of a quincux grid 
// whose edge length is specified by the
// .B -bumpspacing
// option (no default).
// The maximum magnitude of the displacement (in pixels), along any direction is 
// specified by the
// .B -maxdisp
// parameter (no default).
// .I B
// The magnitude of each displacement is the one that maximizes the
// correlation between the displaced 
// .I A
// image 
// and the 
// .I B
// image, in a neighborhood of the bump center.
// .PP
// The two images must have the same width and height.
// The displacement functions will be written to disk as
// two images, called 
// .I name
// .B-dx.
// .Iext
// and 
// .I name
// .B-dy.
// .I ext
// , where 
// .I name 
// is the value of the 
// .B -name
// parameter (defaulting to 
// .B out
// ).
// The extension 
// .I ext
// is 
// .B ppm
// if the input images are in color, and
// .B pgm
// otherwise.
// The output images will have pixels ranging from 
// 0 to 
// .IR depth
// (default 255).
// The displacement values will be multiplied by the 
// .B -scale
// parameter (default 10), and added to the 
// .B -offset
// parameter (default 0), before being rounded to the nearest integer.
// .PP
// All flags can be abbreviated to their shortest unique prefix.
// .SH "SEE ALSO"
// pnmbend(1), pnm(5)
// .SH AUTHOR
// Copyright (C) 2000 by Jorge Stolfi <stolfi@dcc.unicamp.br>.

#include <jspnm.h>
#include <uint16_image.h>
/* #include <sample.h> */
/* #include <jsmisc.h> */
#include <math.h>
#include <stdio.h>

/* Command line arguments */
typedef struct options_t 
  { double bumpspacing;
    double dispmax;
    double scale;
    double offset;
    int depth;
    char *outname;
    char *imAname;
    char *imBname;
    /* Deduced from image types: */
    int chns;  /* 1 or 3 */ 
  } options_t; 

typedef struct ShiftVector
  { double x[3];
    double y[3]; 
  } ShiftVector;  
  /* Shifts per channel and axis at one bump center. */

typedef struct BumpCenter
  { double col;
    double row; 
  } BumpCenter;  
  /* Pixel indices (fractional) of bump center. */

int main(int argc, char **argv);
options_t *get_options(int argc, char **argv);

int main(int argc, char **argv)
  {
    /* Max sampling radius for correlation: */
    double sradius;

    /* Input images: */
    float_image_buffer_t imA; FILE* fpA = NULL;
    float_image_buffer_t imB; FILE* fpB = NULL;

    /* Output images */
    float_image_buffer_t imDX; FILE* fpDX = NULL;
    float_image_buffer_t imDY; FILE* fpDY = NULL;

    BumpCenter  **ctr;     /* Array of bump centers. */
    ShiftVector **shift;   /* Array of shift vectors. */

    /* Parse command line arguments: */
    options_t *o = get_options(argc, argv);

    sradius = o->bumpspacing/2.0;

    open_input_image(o->imAname, &imA, &fpA, sradius);
    open_input_image(o->imBname, &imB, &fpB, sradius + o->maxdisp);
    affirm(imA.cols == imB.cols, "mismatched input image cols");
    affirm(imA.rows == imB.rows, "mismatched input image rows");
    affirm(imA.format == imB.format, "mismatched input image formats");

    compute_shifts(&imA, &imB, sradius, &ctr, &shift);

    open_output_image(cat(cat(o->outname, "-dx"), ext), &imDX, &fpDX);
    open_output_image(cat(cat(o->outname, "-dy"), ext), &imDY, &fpDY);

    compute_shift_images(ctr, shift, &imDX, fpDX, &imDY, fpDY);

    pm_close(fpDX); 
    pm_close(fpDY);
    pm_close(fpA);
    pm_close(fpB);
    return(0);
  }  

void open_input_image(char *name, float_image_buffer_t *im, FILE **fp, double dymax)
  /* Opens file named {name}, saves it in {*fp}.
    Also reads the image header from that file and initializes {*im}.
    Allocates enough buffer rows for sampling the image at pixels {y+dy} where 
    {dy} ranges in {[-dymax .. +dymax]}. */
  {
    /* Compute upper and lower bounds for the displacement DY: */
    int DYmin = floor(-dymax) - FILTER_RADIUS;
    int DYmax = ceil(+dymax) + FILTER_RADIUS;

    fprintf(stderr, "file %s:\n", name);
    fprintf(stderr, "  DYmin = %d DYmax = %d\n", DYmin, DYmax);
    (*fp) = pm_openr(name);
    float_image_read_pnm_buffer_new(*fp, im, DYmin, DYmax, 0.0, 1.0);
    fprintf(stderr, "  rows = %d cols = %d maxval = %d format = %d\n",
      im->rows, im->cols, im->maxval, im->format
    );
  }

void compute_shifts
  ( float_image_buffer_t *imA, float_image_buffer_t *imB,
    double sradius, 
    BumpCenter ***ctrp,
    ShiftVector ***shiftp
  )
  {
    /* Parameters of bump grid: */
    double bDrow = o->bumpspacing * sqrt(3.0)/2.0;
    int brows = floor(((double) imA->rows) / bDrow);
    double bFrow = (((double) imA->rows) - ((double) (brows-1))*bdrow)/2.0;

    double bDcol = o->bumpspacing;
    int bcols = floor(((double) imA->cols) / bDcol);
    double bFcol = (((double) imA->cols) - ((double) (bcols-1))*bDcol)/2.0;

    int nextArow = 0;  /* Next row to be read into image A. */
    int nextBrow = 0;  /* Next row to be read into image B. */

    BumpCenter **ctr = (BumpCenter **)malloc(brows*sizeof(BumpCenter*)); 
    ShiftVector **shift = (ShiftVector **)malloc(brows*sizeof(ShiftVector*));

    affirm (ctr != NULL, "no memory for bump centers");
    affirm (shift != NULL, "no memory for bump shifts");

    for (brow = 0; brow < brows; ++brow)
      { 
        /* Pixel row (fractional) corresponding to bump row */
        double bCrow = bFrow + ((double) brow) * bDrow;

        /* fprintf(stderr, {computing displacements for bump row %d...\n}, brow); */

        /* Read input image lines: */
        ensure_rows_present(imA, fpA, &nextArow, bCrow, sradius);
        ensure_rows_present(imB, fpB, &nextBrow, bCrow, sradius + o->maxdisp);

        /* Allocate bump centers and shift vectors for this bump row: */
        ctr[brow] = (BumpCenter *)malloc(bcols*sizeof(BumpCenter));
        affirm (ctr[brow] != NULL, "no memory for bump centers");
        shift[brow] = (ShiftVector *)malloc(bcols*sizeof(ShiftVector));
        affirm (shift[brow] != NULL, "no memory for bump shifts");

        /* Compute displacements at each bump center: */
        for ( bcol = 0; bcol < bcols; ++bcol)
          {
            float skew = ((double) (brow % 2))/2.0;
            float bCcol = bFcol + (((double) bcol) + skew) * bDcol;
            ctr[brow][bcol] = (BumpCenter){ bCrow, bCcol };
            shiftX[brow][bcol] = compute_shift(imA, imB, sradius, bCrow, bCcol);
          }

      }

    }

  }
  
void ensure_rows_exist(
  float_image_buffer_t *im, File *fp, int *nextRow, double ctrRow, double sradius
)
{
  int minRow = floor(ctrRow - sradius) - FILTER_RADIUS;
  int maxRow = ceil(ctrRow + sradius) + FILTER_RADIUS;
  while ((*nextRow) <= maxRow)
    { float_image_read_pnm_buffer_load_next_row(fp, im, (*nextRow)); (*nextRow)++; }
  affirm(im->inirow <= minRow, "inirow inconsistent");
  affirm(im->finrow >= maxRow, "finrow inconsistent");
}

ShiftVector compute_shift(
  float_image_buffer_t imA, imageBuffer imB, 
  double sradius,
  double bCrow, double BcCol
)
{
  int nch = imA->chns;
  affirm(imA->chns == imB->nChans, "inconsistent chns");
  for( ch = 0; ch <= nch; ch++ )
    { 
      /* Compute shift for channel {ch} */
      ???
    }
}

--- TO FIX ----------------------------------------------------------------------

options_t *parse_args(int argc, char **argv)
  /* Stores in {op} the 
  {

    o->bumpspacing = 0.0;
    o->dispmax = 0.0;
    o->scale = 10.0;
    o->offset = 0.0;
    o->depth = 255;
    o->outname = NULL;
    o->imAname = NULL;
    o->imBname = NULL;


  int argn;


  if ((argn >= argc) || (argv[argn][0] != '-') || (argv[argn][1] == '\0'))
    { arg_error("you must specify \"-bumpspacing\" and \"-maxdisp\""); }
  else if (argparser_keyword_present("-bumpspacing"))
    { argn++;
      parse_shift_args(argc, argv, &argn, &(shiftX[0]), &(shiftY[0]), &shiftchans);
      maptype = MT_FIXED;
    }
  else if (argparser_keyword_present("-maxdisp"))
    { argn++;
      parse_map_args(argc, argv, &argn, &mapDXname, &mapDYname, &offset, &scale);
      maptype = MT_IMAGE;
    }
  else
    { arg_error("bad displacement type"); }


  if (argn == argc) 
    { fpI = stdin; }
  else
    { fpI = pm_openr(argv[argn]); argn++; }

  if (argn != argc) { arg_error("extra arguments"); }
}

 
  /* Open and initialize output images */
  ???
  double *pA;
  double *pB;
  double *pDX;
  double *pDY;
  
  int brow, bcol;        /* Bump row and column indices. */
  double frow, fcol;     /* Fractional pixel row and column of bump center. */
  int row, col;          /* Pixel row and column indices. */
  int rdrow;
  int row, col, r;

    }
    
 
  
  /* Intialize output displacement buffers: */
  float_image_read_pnm_buffer_new(fpI, &imI, DYmin, DYmax, 0.0, 1.0);
  
  /* Check whether displacement maps are consistent with input: */
  if (maptype == MT_IMAGE)
    { int tDX = PNM_FORMAT_TYPE(imDX.format);
      int tDY = PNM_FORMAT_TYPE(imDY.format);
      if (PNM_FORMAT_TYPE(imI.format) != PPM_TYPE)
        { affirm(tDX != PPM_TYPE, "incompatible DX map format");
          affirm(tDY != PPM_TYPE, "incompatible DY map format");
        }
      affirm(imDY.rows == imI.rows, "incompatible displacement map rows");
      affirm(imDY.cols == imI.cols, "incompatible displacement map cols");
    }  
  else
    { if ((PNM_FORMAT_TYPE(imI.format) != PPM_TYPE) && (shiftchans > 1))
        { pm_error("too many shift vectors"); }
    }
  
  /* Initialize output image file: */
  float_image_write_pnm_buffer_new(stdout, &imO, imI.rows, imI.cols, imI.maxval, imI.format, 0, 0, 0.0, 1.0);

  fprintf(stderr, "opened output image  rows = %d cols = %d maxval = %d format = %d\n",
    imO.rows, imO.cols, imO.maxval, imO.format
  );
  
}

      float_image_write_pnm_buffer_dump_first_row(stdout, &imO, row);
      /* fprintf(stderr, "done.\n"); */

--- JUNK ----------------------------------------------------------------------

extern double strtod ARGS((const char *, char **));
void arg_error ARGS((char *msg));
void parse_shift_args ARGS((int argc, char **argv, int *argnp, double shiftX[], double shiftY[], int *nchansp));
void parse_map_args ARGS((int argc, char **argv, int *argnp, char **mapDXnamep, char **mapDYnamep, double *offsetp, double *scalep));

typedef enum {MT_FIXED = 0, MT_IMAGE = 1} MapType;


    argparser_t pp = argparser_new(stderr, argc, argv);

    argn = 1;


void parse_shift_args(int argc, char **argv, int *argnp, double *shiftX, double *shiftY, int *nchansp)
{ 
  /* Parses the arguments for the "-shift" mode. */
  /* Sets shiftX, shiftY. */
  
  int argn = (*argnp);
  int r; 
  if (argn+6 < argc)
    { char *restX, *restY;
      /* fprintf (stderr, "three-chanel shift (%s .. %s)\n", argv[argn], argv[argn+5]); */
      (*nchansp) = 3;
      for (r = 0; r < 3; r++)
        { shiftX[r] = strtod(argv[argn++], &restX);
          shiftY[r] = strtod(argv[argn++], &restY);
          if ((*restX != '\0') || (*restY != '\0')) { arg_error("bad shift"); } 
        }
    }
  else if (argn+2 < argc)
    { char *restX, *restY;
      /* fprintf (stderr, "single-chanel shift (%s %s)\n", argv[argn], argv[argn+1]); */
      (*nchansp) = 1;
      shiftX[0] = strtod(argv[argn++], &restX); 
      shiftX[1] = shiftX[0]; shiftX[2] = shiftX[0];
      shiftY[0] = strtod(argv[argn++], &restY);
      shiftY[1] = shiftY[0]; shiftY[2] = shiftY[0];
      if ((*restX != '\0') || (*restY != '\0')) { arg_error("bad shift"); } 
    }
  else 
    { arg_error("not enough arguments for \"-shift\""); }
  (*argnp) = argn;
}

void parse_map_args(int argc, char **argv, int *argnp, char **mapDXnamep, char **mapDYnamep, double *offsetp, double *scalep)
{ 
  /* Parses the arguments for the "-map" mode. */
  /* Sets mapXname, mapYname, scale, offset. */
  
  int argn = (*argnp);

  if (argn+2 < argc)
    { (*mapDXnamep) = argv[argn+1]; 
      (*mapDYnamep) = argv[argn+2]; 
      argn += 2;
    }
  else
    { arg_error("not enough arguments for \"-map\""); }

  (*offsetp) = 0.0;
  (*scalep) = 1.0;
  while ((argn < argc) && (argv[argn][0] == '-') && (argv[argn][1] != '\0'))
    { if (argparser_keyword_present("-scale"))
        { if (argn+1 < argc)
            { char *rest;
              (*scalep) = argparser_get_next_double(pp, ???, ???); ++argn;
              if (*rest != '\0') { arg_error("bad scale factor"); } 
            }
          else
            { arg_error("missing scale factor"); }
        }
      else if (argparser_keyword_present("-offset"))
        { if (argn+1 < argc)
            { 
              char *rest;
              (*offsetp) = argparser_get_next_double(pp, ???, ???); ++argn;
              if (*rest != '\0') { arg_error("bad offset"); } 
            }
          else
            { arg_error("missing offset"); }
        }
      else
        { arg_error("unrecognized option"); }
      ++argn;
    }
  (*argnp) = argn;
}
  
void arg_error(char *msg)
{
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
  
    PROG_NAME " version " PROG_VERS ", usage:\n"
    PROG_HELP; 
  pm_message(msg);
  pm_usage(usage);
}

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 2000 by the State University of Campinas (UNICAMP).
**
** Created sep/2000 by Jorge Stolfi, IC-UNICAMP.       
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
