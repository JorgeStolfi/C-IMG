#define PROG_NAME "fni_std_extract"
#define PROG_DESC "converts FNI images to a standard size and shape"
#define PROG_VERS "1.0"

/* !!! SHOULD BE MERGED INTO {plan_extract.c} !!! */

/* Last edited on 2021-12-31 23:42:35 by stolfi */

#define fni_std_extract_C_COPYRIGHT \
  "Copyright © 2007 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -minSize {EX0} {EY0} \\\n" \
  "  -midSize {EX1} {EY1} \\\n" \
  "  -maxSize {EXM} {EYM} \\\n" \
  "  -maxScale {MAXSCALE} \\\n" \
  "  -maxStretch {MAXSTRETCH} \\\n" \
  "  -maxPad {MAXPAD} \\\n" \
  "  -maxCut {MAXCUT} \\\n" \
  "  -prefix {PREFIX} \\\n" \
  "  < {INPUT_FNI}"

#define PROG_WHAT \
  "  Reads a FNI file from standard input, and extracts" \
  " from it one or more sub-images, by trimming and/or scaling" \
  " and/or cropping and/or padding.\n" \
  "\n" \
  "  For this description, the {X} axis is horizontal from the" \
  " left edge, and the {Y} axis is vertical from the bottom" \
  " edge.\n" \
  "\n" \
  "  In the first stage, any background margins are identified" \
  " and trimmed.  More precisely, each edge of the input" \
  " image {I} is examined to see whether its pixels have" \
  " uniform color {C} (apart from small levels of noise" \
  " or a few exceptions).  If so, that edge is marked" \
  " as \"paddable\" with color {C}.  The program then" \
  " removes any rows or columns of pixels along" \
  " paddable edges of {I} that are uniformly colored" \
  " with the pad color of that edge.\n" \
  "\n" \
  "  Let {O} be the image reslting from this trimming," \
  " and {OX × OY} its dimensions.  The image {O} is then" \
  " scaled down to an image {S} with some size" \
  " {SX × SY}.  Then {S} is cropped and/or padded to produce" \
  " one or more images {E[r,s]}, for {r,s} ranging over some" \
  " integer intervals {rmin..rmax} and {smin..smax} that straddle 0.\n" \
  "\n" \
  "  All extracted images {E[r,s]} have the same dimensions" \
  " {EX × EY}.  Image {E[0,0]} is concentric with {S} (rounding" \
  " down if either has odd size), and image {E[r,s]} is" \
  " displaced from {E[0,0]} by {r*EX/2} in the {X} direction" \
  " and {s*EY/2} in the {Y} direction.\n" \
  "\n" \
  "  The extracted image dimensions {EX × EY} belong to a finite" \
  " geometric progression with ratio 2, or to two such" \
  " progressions, interleaved, which are specified through" \
  " the options \"-minSize\", \"-midSize\", and \"-maxSize\".\n" \
  "\n" \
  "  For each extracted image {E[r,s]}, the program" \
  " writes out a FNI file called \"{PFX}-{R}{S}.fni\", where {PFX} is" \
  " the \"-profix\" parameter, {R} is the digit {r+5}," \
  " and {S} is the digit {s+5}. (The program fails if {|r|}" \
  " or {|s|} exceed 4).\n" \
  "\n" \
  "  For each extracted image {E[r,s]}, the program also writes a" \
  " text file called \"{PFX}-{R}{S}.txt\", that" \
  " describes how the image was extracted.  The file" \
  " has the format\n" \
  "\n" \
  "     size {IX}x{IY} cropped to {OX}x{OY}+{IODX}+{IODY}" \
  " scaled to {SX}x{SY} segment {EX}x{EY}±{SEDX}±{SEDY}\n" \
  "\n" \
  "  where {IODX,IODY} is the displacements of {O} relative" \
  " to {I}, and {SEDX,SEDY} is the displacement of {E[r,s]} relative" \
  " to {S}.  Note that displacements here refer to the UPPER LEFT" \
  " corners of the two images, with the {X} axis pointing right" \
  " but the {Y} axis pointing DOWN.  If the second image" \
  " is not contained in the first, it is understood that" \
  " the overflowing edges are filled with the padding (dominant)" \
  " color of the corresponding edge of {I}.\n"
  
#define PROG_OPTS \
  "  -minSize {EX0} {EY0}\n" \
  "    Specifies the minimum dimensions of an extracted" \
  " image.  This parameter is mandatory.\n" \
  "\n" \
  "  -midSize {EX1} {EY1}\n" \
  "    Specifies the second largest allowed dimensions" \
  " for an extracted image.  The dimensions must lie between {EX0 × EY0} and" \
  " {2*EX0 × 2*EY0}, and the ratio {EX1/EX0} must be equal to {EY1/EY0}.  If this" \
  " parameter is not specified, or if it is equal to" \
  " either {EX0 × EY0} or {2*EX0 × 2*EY0}, the program" \
  " assumes that there are no intermediate image dimensions.\n" \
  "\n" \
  "  -maxSize {EXM} {EYM}\n" \
  "    Specifies the maximum allowed size of an" \
  " extracted image.  This parameter is mandatory.\n" \
  "\n" \
  "  -rotate\n" \
  "    If present, this option specifies that the shape" \
  " of the extracted images be turned 90 degrees.  In other" \
  " words, if the extracted images are allowed to have" \
  " dimensions {EX × EY}, they are allowed to have dimensions" \
  " {EY × EX} too.  If this option is not specified," \
  " the extracted images must have the specified shape," \
  " including orientation.\n" \
  "\n" \
  "  -maxScale {MAXSCALE}\n" \
  "    Specifies the maximum value for the scaling factor {SX/OX = SY/OY}" \
  " when converting from image {O} to image {S}.\n" \
  "\n" \
  "  -maxStretch {MAXSTRETCH}\n" \
  "    Specifies the maximum relative stretching factor of the" \
  " rectangle {SX × SY}, relative to a rectangle with the same" \
  " proportions as {O}, along its longest side.  The rectangle" \
  " {SX × SY} may also be shrunk by {1/MAXSTRETCH} along that" \
  " side.\n" \
  "\n" \
  "  -maxPad {MAXPAD}\n" \
  "    Specifies the maximum fraction of an extracted image" \
  " that may be filled by padding, along either axis.  If" \
  " omitted, {MAXPAD} is assumed to be 0.\n" \
  "\n" \
  "  -maxCut {MAXCUT}\n" \
  "    Specifies the maximum fraction of the scaled" \
  " image {S} that may remain uncovered by the extracted" \
  " images.  If omitted, {MAXCUT} is assumed to be 1.\n" \
  "\n" \
  "  -prefix {PREFIX}\n" \
  "    Specifies the prefix for output file names. This" \
  " parameter is mandatory.\n" \

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_WHAT

#define PROG_FOOO \
  PROG_UGGH \
  "\n" \
  "OPTIONS\n" \
  PROG_OPTS \
  "\n" \
  argparser_help_info_INFO "\n" \
  "\n" \
  "SEE ALSO" \
  "  pnmxpad(1)" \
  "AUTHOR" \
  "  Created on 2007-04-06 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " fni_std_extract_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
 
/* We need to set these in order to get {isnan}. What a crock... */
#undef __STRICT_ANSI__
#define _ISOC99_SOURCE 1
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <bool.h>
#include <i2.h>
#include <irange.h>
#include <argparser.h>

#include <frgb.h>
#include <frgb_ops.h>
#include <float_image.h>
#include <float_image_trim.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t 
  { double maxScale;   /* Maximum {S/O} linear size ratio. */
    double maxStretch; /* Maximum aspect stretch factor in {O->S} scaling. */
    double maxPad;     /* Maximum linear fraction of {E[r,s]} not covered by {S}. */
    double maxCut;     /* Maximum width fraction of {S} covered by {E[r,s]}. */
    /* Parameters that describe the valid sizes for extracted images {E[r,s]}: */
    i2_t minSize;      /* Minimum dimensions in main series. */
    i2_t midSize;      /* Minimum dimesnions in intermediate series. */
    i2_t maxSize;      /* Maximum dimensions for either series. */
    bool_t rotate;     /* TRUE considers also the rotated shapes. */
    /* Parameters for background identification: */
    double noise;      /* Maximum deviation from average color. */
    double except;     /* Maximum fraction of exceptional pixels. */
    /* Output parameters: */
    char *prefix;      /* Prefix for output file names. */
  } options_t;

/* INTERNAL TYPES */
 
/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

options_t *fse_get_options(int argc, char **argv);

void process_image(float_image_t *img, options_t *o);
  /* Applies the analysis to the image {img},
    writes the extracted images {E[r,s]} to files
    "{o.prefix}-{R}{S}.txt" and "{o.prefix}-{R}{S}.fni"
    where {R = 5+r} and {S = 5+s}. */

float_image_t *fse_read_float_image(char *name);
  /* Reads a float image, in the format of {float_image_write}, from
    the named file (or {stdin} if {name} is "-"). */
    
void fse_write_float_image(char *name, float_image_t *img);
  /* Writes the image {img}, in the format of {float_image_write}, to
    the named file (or to {stdout} if {name} is "-"). */

void fse_error(char *msg);
  /* Aborts the program with error message {msg}. */

float_image_t *fse_trim_and_scale(float_image_t *I, irange_t OBox[], i2_t *SSize);
  /* Extracts the sub-image of {I} contained in the rectangle {OBox}, and
    scales it down to size {SSize}. */

void fse_extract_and_write_image(float_image_t *S, int r, int s, i2_t *ESize, i2_t *EPos, char *prefix);
  /* Extracts the sub-image {E[r,s]} from image {S} and writes it out
    to a FNI file called "{prefix}-{R}{S}.fni" where {R} is the digit
    {r+5} and {S} is the digit {s+5}.
    
    The sub-image starts at position {EPos + (r*DX,s*DY)} where
    {DX=floor(ESize.c[0]/2)} and {DY=floor(ESize.c[1]/2)}. The
    displacement refers to the LOWER LEFT corners with Y axis pointing UP.
    
    Fails if {|r|} or {|s|} exceed 4. */
    
i2_t fse_parse_size(argparser_t *pp, int MAX);
  /* Parses the next two arguments from the command line as two integers {NX} {NY},
    each between 1 and {MAX}.  Returns the pair {(NX,NY)}. */

/* Maximum rows, columns, elements (to avoid absurd allocs): */
#define fse_MAX_PIXELS (16*256*256*256)   
#define fse_MAX_SIZE fse_MAX_PIXELS

#define TwoPi (2.0 * M_PI)

int main_debug = TRUE;

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  { 
    /* Parse command line arguments: */
    options_t *o = fse_get_options(argc, argv);

    /* Read image: */
    float_image_t *img = fse_read_float_image("-");
    
    /* Process image: */
    process_image(img, o);

    fprintf(stderr, "done.\n");
    return 0;
  }
  
void process_image(float_image_t *I, options_t *o)
  {
    /* Get raw image dimensions: */
    int NC, IX, IY;
    float_image_get_size(I, &NC, &IX, &IY);
    i2_t ISize = (i2_t){{IX, IY}}; /* Dimensions of {I}. */

    /* First stage: trim margins and remember pad colors, if any. */
    bool_t paddable[4]; /* TRUE if edge {e} is paddable. */
    frgb_t padcolor[4]; /* Pad color along each edge, or {(-1,-1,-1)} if unpaddable. */
    irange_t OBox[2]; /* Domain of image {O}, relative to {I}. */
    float_image_trim_background(I, o->noise, o->except, OBox, paddable, padcolor);
    /* Get dimensions {OX,OY} of {O}: */
    int OX = OBox[0].end[1] - OBox[0].end[0];
    int OY = OBox[1].end[1] - OBox[1].end[0];
    i2_t OSize = (i2_t){{OX, OY}};
    
    /* Find the best choice for {SX,SY,EX,EY}: */
    plan_t pl = fse_find_best_plan(&OSize, paddable, o);
    if (pl.score < 0) { fse_error("no acceptable plan"); }
    
    /* Scale down the input image: */
    float_image_t *S = fse_trim_and_scale(I, OBox, &(pl.SSize));
    
    /* Extract the sub-images and write them out: */
    int rmin = pl.tr[0].end[0], rmax = pl.tr[0].end[1];
    int smin = pl.tr[1].end[0], smax = pl.tr[1].end[1];
    int r, s;
    for (r = rmin; r <= rmax; r++)
      for (s = smin; s <= smax; s++)
        { /* Extract tile {E[r,s]} and write it out: */
          fse_extract_and_write_image(S, r, s, &(pl.ESize), &(pl.EPos), o->prefix);
          fse_write_description(&ISize, OBox, &(pl.SSize), r, s, &(pl.ESize), &(pl.EPos), o->prefix); 
        }
  }
  
options_t *fse_get_options(int argc, char **argv)
  {
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    argparser_get_keyword(pp, "-minSize");
    o->minSize = fse_parse_size(pp, fse_MAX_SIZE); 
    
    if (argparser_keyword_present(pp, "-midSize"))
      { o->midSize = fse_parse_size(pp, fse_MAX_SIZE);
        int EX0 = o->minSize.c[0], EY0 = o->minSize.c[1];
        int EX1 = o->midSize.c[0], EY1 = o->midSize.c[1];
        /* Check whether {o.midSize} is in the proper range: */
        if ((EX1 < EX0) || (EX1 > 2*EX0) || (EY1 < EY0) || (EY1 > 2*EY0))
          { argparser_error(pp, "min intermediate dimensions are out of range"); }
        /* Check whether {o.midSize} has the same aspect as {o.minSize}: */
        if (EX0*EY1 != EY0*EX1)
          { argparser_error(pp, "min intermediate dimesnions have wrong aspect"); }
      } 
    else
      { /* No intermediate series: */
        o->midSize = o->minSize;
      }
      
    argparser_get_keyword(pp, "-maxSize");
    o->maxSize = fse_parse_size(pp, fse_MAX_SIZE); 
    
    if (argparser_keyword_present(pp, "-maxScale"))
      { o->maxScale = argparser_get_next_double(pp, 0.0, 1.0); }
    else 
      { o->maxScale = 1.0; }
      
    if (argparser_keyword_present(pp, "-maxPad"))
      { o->maxPad = argparser_get_next_double(pp, 0.0, 1.0); }
    else
      { o->maxPad = 1.0; }
      
    if (argparser_keyword_present(pp, "-maxCut"))
      { o->maxCut = argparser_get_next_double(pp, 0.0, 1.0); }
    else
      { o->maxCut = 1.0; }
      
    if (argparser_keyword_present(pp, "-noise"))
      { o->noise = argparser_get_next_double(pp, 0.0, 1.0); }
    else
      { o->noise = 0.02; }
      
    if (argparser_keyword_present(pp, "-except"))
      { o->except = argparser_get_next_double(pp, 0.0, 1.0); }
    else
      { o->except = 0.05; }
      
    argparser_get_keyword(pp, "-prefix");
    o->prefix = argparser_get_next(pp); 
    
    /* Check for extraneous arguments: */
    argparser_finish(pp);
        
    return o;
  }
  
i2_t fse_parse_size(argparser_t *pp, int MAX)
  { i2_t sz;
    sz.c[0] = argparser_get_next_int(pp, MAX);
    sz.c[1] = argparser_get_next_int(pp, MAX);
    return sz; 
  }

