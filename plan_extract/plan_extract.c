#define PROG_NAME "plan_extract"
#define PROG_DESC "selects a sub-image with standard size and shape"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-20 06:13:09 by stolfi */

#define plan_extract_C_COPYRIGHT \
  "Copyright © 2007 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -input {PREFIX} {IX} {IY} \\\n" \
  "  [ -take {TX} {TY} {ITDX} {ITDY} ] \\\n" \
  "  [ -pad {SIDES} ] \\\n" \
  "  -minSize {EX0} {EY0} \\\n" \
  "  [ -midSize {EX1} {EY1} ] \\\n" \
  "  -maxSize {EXM} {EYM} \\\n" \
  "  [ -rotate ] \\\n" \
  "  [ -maxScale {MAXSCALE} ] \\\n" \
  "  [ -maxStretch {MAXSTRETCH} ] \\\n" \
  "  [ -maxPad {MAXPAD} ] \\\n" \
  "  [ -maxCut {MAXCUT} ] \\\n" \
  "  [ -debug {LEVEL} ]"

#define PROG_WHAT \
  "  Selects one or more sub-images from some image, by" \
  " cropping and/or padding and/or scaling.  The" \
  " program does not do the image extraction" \
  " itself; it only computes the parameters for standard" \
  " image cropping/padding/scaling tools.\n" \
  "\n" \
  "  For this description, the {X} axis is horizontal from the" \
  " left edge, and the {Y} axis is vertical from the bottom" \
  " edge.\n"

#define PROG_STGS \
  "  The program assumes that the input image {I} has given" \
  " dimensions {IX} by {IY}, and that the sub-image(s) will be" \
  " extracted from it in four stages:\n" \
  "   * A sub-image {C} is cropped from {I};" \
  "   * An image {P} is obtained from {C} by padding along some edges;" \
  "   * An image {S} is obtained from {P} by scaling." \
  "   * One or more images {E[r,s]} are cropped from {S}." \
  "\n" \
  "  All extracted images {E[r,s]} have the same dimensions" \
  " {EX × EY}, selected by the program, and fit inside {S}." \
  "  The indices {r,s} will range over some integer intervals" \
  " {rmin..rmax} and {smin..smax} that straddle 0.  Image" \
  " {E[0,0]} will be concentric with {S}, and image {E[r,s]} will be" \
  " displaced from {E[0,0]} by {r*EX/2} in the {X} direction" \
  " and {s*EY/2} in the {Y} direction.  Thus, extracted" \
  " images with successive indices will have 50% overlap\n" \
  "\n" \
  "  For each extracted image {E[r,s]}, the program writes a" \
  " text file called \"out/{PREFIX}-{RR}-{SS}.txt\", that" \
  " describes how the image was extracted; where {RR} and {SS} are" \
  " two-digit values of {r+50} and {s+50}.  The file" \
  " has the format\n" \
  "\n" \
  "     size {IX} {IY}\n" \
  "     cropped to {CX} {CY} {ICDX} {ICDY}\n" \
  "     padded to {PX} {PY} {CPDX} {CPDY}\n" \
  "     scaled to {SX} {SY}\n" \
  "     segment {EX} {EY} {SEDX} {SEDY}\n" \
  "     score {SCORE}\n" \
  "\n" \
  " where\n" \
  "\n" \
  "     {CX} {CY} {ICDX} {ICDY}  is the domain of {C} relative to {I};\n" \
  "     {PX} {PY} {CPDX} {CPDY}  is the domain of {P} relative to {C};\n" \
  "     {SX} {SY}                is the size of {S};\n" \
  "     {EX} {EY} {SEDX} {SEDY}  is the domain of {E[r,s]} relative to {S}; and\n" \
  "     {SCORE}                  is a quality score for this plan (1 best, 0 worst).\n" \
  "\n" \
  " Note that all displacements here refer to the UPPER LEFT" \
  " corners of the two images, with the {X} axis pointing right" \
  " but the {Y} axis pointing DOWN.  If the second image" \
  " is not contained in the first, it is understood that" \
  " the overflow areas are to be padded with some appropriate color.\n"
  
#define PROG_PSEL \
  "  The parameters of all three stages are chosen by the program, based" \
  " on preferences and constraints specified on the command" \
  " line.  They are\n" \
  "   * the dimensions {CX × CY} of the cropped image {C};\n" \
  "   * the coordinates {(ICDX,ICDY)} of the upper left corner of {C}, relative to {I};\n" \
  "   * the dimensions {PX × PY} of the padded image {P};\n" \
  "   * the coordinates {(CPDX,CPDY)} of the lower left corner of {P}, relative to {C};\n" \
  "   * the dimension {SX × SY} of the scaled image {S}; and\n" \
  "   * the dimensions {EX × EY} of each output image.\n" \
  "\n" \
  "  The extracted image dimensions {EX × EY} belong to a finite" \
  " geometric progression with ratio 2, or to two such" \
  " progressions, interleaved, which are specified through" \
  " the options \"-minSize\", \"-midSize\", and \"-maxSize\".\n" \
  "\n" \
  "  The cropped image {C} is contained in {I} (or in a" \
  " sub-image {T} of {I} specified by the user).\n"
  
#define PROG_OPTA \
  "  -input {PREFIX} {IX} {IY}\n" \
  "    Specifies that the program should consider a single" \
  " image {I} with {IX} columns and {IY} rows.  The {PREFIX} is" \
  " used to form the output file names.  This option" \
  " is mandatory.\n" \
  "  -pad {SIDES}\n" \
  "    Specifies the sides along which the sub-image {T}" \
  " can be padded (with some unspecified background color).  The" \
  " {SIDES} parameter is a string that consists" \
  " of any combination of the letters 'l' (left), 'r' (right), 't' (top)," \
  " or 'b' (bottom), in any order.  If this option is omitted," \
  " the program assumes \"-pad tblr\" (i.e. {T} can be padded" \
  " on all edges).\n" \
  "\n" \
  "  -take {TX} {TY} {ITDX} {ITDY}\n" \
  "    Specifies the valid sub-rectangle of the input" \
  " image.  The rectangle has {TX} columns and {TY}" \
  " rows; it starts {ITDX} columns to the right of the" \
  " left edge of {I} and {ITDY} rows below the top edge" \
  " of {I}.  If this parameter is not specified, the program" \
  " assumes that the whole image is valid," \
  " namely \"-take {IX} {IY} 0 0\".\n" \
  "\n" \
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
  " of the extracted images may be turned 90 degrees.  In other" \
  " words, if the extracted images are allowed to have" \
  " dimensions {EX × EY}, they are allowed to have dimensions" \
  " {EY × EX} too.  If this option is not specified," \
  " the extracted images must have the specified shape," \
  " including orientation.\n" \

#define PROG_OPTB \
  "  -maxScale {MAXSCALE}\n" \
  "    Specifies the maximum value for the scaling factor {SX/CX = SY/CY}" \
  " when converting from image {P} to image {S}.  It must be positive" \
  " and no greater than 1.  If this option is omitted," \
  " {MAXSCALE} is assumed to be 1 (i.e. the image may be" \
  " reduced but not enlarged).\n" \
  "\n" \
  "  -maxStretch {MAXSTRETCH}\n" \
  "    Specifies the maximum relative stretching factor of the" \
  " rectangle {SX × SY}, relative to a rectangle with the same" \
  " proportions as {O}, along its longest side.  The rectangle" \
  " {SX × SY} may also be shrunk by {1/MAXSTRETCH} along that" \
  " side.  If this option is omitted," \
  " {MAXSTRETCH} is assumed to be 0 (i.e. the {P}" \
  " and {S} images must have precisely the same aspect ratio).\n" \
  "\n" \
  "  -maxPad {MAXPAD}\n" \
  "    Specifies the maximum fraction of any extracted image" \
  " that may taken up by such padding, along either" \
  " axis.  If this option is omitted, {MAXPAD} is assumed to" \
  " be 0 (i.e. padding is not allowed).\n" \
  "\n" \
  "  -maxCut {MAXCUT}\n" \
  "    Specifies the maximum fraction of the scaled" \
  " image {S} that may remain uncovered by the extracted" \
  " images.  If this option is omitted, {MAXCUT} is" \
  " assumed to be 0.5 (i.e. up to 50% of the {T} image" \
  " may be discarded).\n" \
  "\n" \
  "  -debug {LEVEL}\n" \
  "    Specifies the verbosity of debugging" \
  " printouts.  The {LEVEL} is an integer. If {LEVEL <= 0}, only error messages" \
  " are printed.  Typically each increase by 10 in {LEVEL} generates" \
  " debugging output for one additional level of procedure call.\n" \

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_WHAT \
  "\n" \
  "SUB-IMAGE EXTRACTION PROCESS\n" \
  PROG_STGS \
  "\n" \
  "PARAMETER SELECTION\n" \
  PROG_PSEL \
  "\n" \
  "OPTIONS\n" \
  PROG_OPTA \
  "\n" \
  PROG_OPTB \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO" \
  "  pnmxpad(1)" \
  "AUTHOR" \
  "  Created on 2007-05-05 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " plan_extract_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
 
/* We need to set these in order to get {isnan}. What a crock... */
#undef __STRICT_ANSI__
#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <jsfile.h>
#include <i2.h>
#include <irange.h>
#include <argparser.h>

#include <ple_find.h>

/* COMMAND-LINE OPTIONS */
  
typedef struct options_t 
  { /* Input image: */
    char *prefix;        /* Prefix for output file names. */
    i2_t ISize;          /* Input image size. */
    /* Cropping parameters: */
    irange_t TBox[2];    /* Limiting bounding box for cropping. */
    /* Padding parameters: */
    bool_t paddable[4];  /* Indicates which edges of {TBOX} may be padded. */
    /* Constraints: */
    double maxScale;     /* Maximum {S/C} linear size ratio. */
    double maxStretch;   /* Maximum aspect stretch factor in {C->S} scaling. */
    double maxPad;       /* Maximum linear fraction of {E[r,s]} not covered by {S}. */
    double maxCut;       /* Maximum width fraction of {S} covered by {E[r,s]}. */
    /* Parameters that describe the valid sizes for extracted images {E[r,s]}: */
    i2_t minSize;        /* Minimum dimensions in main series. */
    i2_t midSize;        /* Minimum dimesnions in intermediate series. */
    i2_t maxSize;        /* Maximum dimensions for either series. */
    bool_t rotate;       /* TRUE considers also the rotated shapes. */
    /* Debugging parameters: */
    int debug;           /* Debugging level: positive means print debug info. */
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

options_t *ple_get_options(int argc, char **argv);

plan_t ple_find_best_plan(options_t *o);
  /* Finds the best image extraction plan, given the command line options {o}. */

void ple_error(char *msg);
  /* Aborts the program with error message {msg}. */

void ple_write_plan(i2_t *ISize, plan_t *pl, int r, int s, char *prefix);
  /* Writes the description of image {E[r,s]}, extracted according to
    plan {pl}, to a text file, in the format explained in the documentation.
    
    Assumes the sub-image {E[r,s]} starts at position {pl->EPos +
    (r*pl.EStep.c[0], s*pl.EStep.c[1])} relative to {S}.
    
    When calling this procedure, the displacements in {pl} refer to
    the LOWER LEFT corners with Y axis pointing UP. In the output
    file, the displacements refer to the UPPER LEFT corners, with the
    Y axis pointing DOWN.
    
    The text file will be called "out/{prefix}-{RR}-{SS}.txt" where {RR}
    is the two-digit integer {r+50} and {SS} is the two-digit integer
    {s+50}. Fails if {|r|} or {|s|} exceed 49. */

i2_t ple_parse_size(argparser_t *pp, int MAX);
  /* Parses the next two arguments from the command line as two integers {NX} {NY},
    each between 1 and {MAX}.  Returns the pair {(NX,NY)}. */

i2_t ple_parse_pos(argparser_t *pp, int MAX);
  /* Parses the next two arguments from the command line as two integers {DX} {DY},
    each between {-MAX} and {+MAX}.  Returns the pair {(DX,DY)}. */

void ple_parse_edge_set(argparser_t *pp, bool_t es[]);
  /* Parses the next argument from the command line as a set of edges,
    where each edge is indicated by a letter 'l', 'r', 't', or 'b',
    in any order. */

void ple_parse_edge_set(argparser_t *pp, bool_t es[]);

int main_debug = TRUE;

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  { 
    /* Parse command line arguments: */
    options_t *o = ple_get_options(argc, argv);

    /* Find the best image-extraction plan: */
    plan_t pl = ple_find_best_plan(o);
    if (pl.score < 0) { ple_error("no acceptable plan"); }
    
    /* Write the extraction data: */
    int rmin = pl.tr[0].end[0], rmax = pl.tr[0].end[1];
    int smin = pl.tr[1].end[0], smax = pl.tr[1].end[1];
    int r, s;
    for (r = rmin; r <= rmax; r++)
      for (s = smin; s <= smax; s++)
        { ple_write_plan(&(o->ISize), &pl, r, s, o->prefix); }
    
    fprintf(stderr, "done.\n");
    return 0;
  }
  
plan_t ple_find_best_plan(options_t *o)
  {
    plan_t pl_best;
    
    bool_t verbose = (o->debug > 0);
    if (verbose) 
      { fprintf(stderr, "+ ple_find_best_plan\n"); 
        ple_debug_box(stderr, "  TBox = ", o->TBox, "\n");
        fprintf(stderr, "\n");
      }
    
    
    pl_best.score = -INF;
    bool_t mid;  /* FALSE for main size series, TRUE for intermediate series. */
    bool_t turn; /* FALSE for unturned sizes, TRUE for turned sizes. */
    for (mid = FALSE; mid <= TRUE; mid++)
      for (turn = FALSE; turn <= o->rotate; turn++)
        { i2_t minESize = (mid ? o->midSize : o->minSize);
          i2_t maxESize = o->maxSize;
          if (turn) 
            { /* Swap axes in {minESize}, {fin{Size}: */
              int t;
              t = minESize.c[0]; minESize.c[0] = minESize.c[1]; minESize.c[1] = t;
              t = maxESize.c[0]; maxESize.c[0] = maxESize.c[1]; maxESize.c[1] = t;
            }
          /* Enumerate valid image sizes {EX,EY} and consider each in turn: */
          i2_t ESize = minESize;
          while ((ESize.c[0] <= maxESize.c[0]) && (ESize.c[1] <= maxESize.c[1]))
            { /* The displacement between tiles is {ESize/2}, rounded down: */
              i2_t EStep = (i2_t){{ ESize.c[0]/2, ESize.c[1]/2 }};
              /* Consider plans that extract images of size {ESize}: */
              ple_eval_plans_for_esize
                ( o->TBox, o->paddable, &ESize, &EStep,
                  o->maxScale, o->maxStretch, o->maxPad, o->maxCut, 
                  &pl_best, o->debug - 10
                );
              /* Get next size pair in same series: */
              ESize.c[0] *= 2; ESize.c[1] *= 2;
            }
        }
    if (verbose) 
      { ple_debug_plan(stderr, "best plan:\n", o->TBox, &pl_best, "\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "- ple_find_best_plan\n");
      }
    return pl_best;
  }


void ple_error(char *msg)
  { fprintf(stderr, "** %s\n", msg);
    exit(1);
  }

void ple_write_plan(i2_t *ISize, plan_t *pl, int r, int s, char *prefix)
  {
    demand(abs(r) <= 49, "bad {r}");
    demand(abs(s) <= 49, "bad {s}");
    
    /* Construct file name: */
    char *fname = jsprintf("out/%s-%02d-%02d.txt", prefix, r+50, s+50);
    FILE *wr = open_write(fname, TRUE);
    
    /* Sizes of {C} and {P} sub-images: */
    i2_t CSize = ple_box_size(pl->CBox);
    i2_t PSize = ple_box_size(pl->PBox);

    /* Print size of initial image {I}: */
    int IX = ISize->c[0], IY = ISize->c[1];
    fprintf(wr, "size %d %d\n", IX, IY);
    
    /* Print domain of crop image {C} relative to {I}: */
    int CX = CSize.c[0], CY = CSize.c[1];
    int ICDX = pl->CBox[0].end[0], ICDY = IY - pl->CBox[1].end[1];
    fprintf(wr, "cropped to %d %d %+d %+d\n", CX, CY, ICDX, ICDY);
    
    /* Print domain of pad image {P} relative to {C}: */
    int PX = PSize.c[0], PY = PSize.c[1];
    int CPDX = pl->PBox[0].end[0], CPDY = CY - pl->PBox[1].end[1];
    fprintf(wr, "padded to %d %d %+d %+d\n", PX, PY, CPDX, CPDY);
    
    /* Print size of scaled image {S}: */
    int SX = pl->SSize.c[0], SY = pl->SSize.c[1];
    fprintf(wr, "scaled to %d %d\n", SX, SY);
    
    /* Print domain of segment {E[r,s]} relative to {S}: */
    int EX = pl->ESize.c[0], EY = pl->ESize.c[1];
    int SEDX = pl->EPos.c[0] + r*pl->EStep.c[0];
    int SEDY = SY - (pl->EPos.c[1] + pl->ESize.c[1] + s*pl->EStep.c[1]);
    fprintf(wr, "segment %d %d %+d %+d\n", EX, EY, SEDX, SEDY);
    
    /* Print score of segment {E[r,s]}: */
    double segScore = pl->score;
    if (segScore > 0) { segScore *= pow(1/M_SQRT2, abs(r) + abs(s)); }
    fprintf(wr, "score %8.6f\n", segScore);

    fclose(wr);
    free(fname);
  }

options_t *ple_get_options(int argc, char **argv)
  {
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    int ax;
    
    argparser_get_keyword(pp, "-input");
    o->prefix = argparser_get_next(pp); 
    o->ISize = ple_parse_size(pp, ple_MAX_SIZE); 
    
    if (argparser_keyword_present(pp, "-take"))
      { /* Parse the relevant box size {TSize}: */
        i2_t TSize = ple_parse_size(pp, ple_MAX_SIZE); 
        /* Parse graphics-style relevant box position {TPos} (upper left, {Y} down): */
        i2_t TPos = ple_parse_pos(pp, ple_MAX_SIZE); 
        /* Convert to math-style position (lower left, {Y} up): */
        o->TBox[0].end[0] = TPos.c[0];
        o->TBox[0].end[1] = o->TBox[0].end[0] + TSize.c[0];
        o->TBox[1].end[1] = o->ISize.c[1] - TPos.c[1];
        o->TBox[1].end[0] = o->TBox[1].end[1] - TSize.c[1];
        /* Check whether {TBox} lies inside {I}: */
        for (ax = 0; ax < 2; ax++)
          { if ((o->TBox[ax].end[0] < 0) || (o->TBox[ax].end[1] > o->ISize.c[ax])) 
              { argparser_error(pp, "\"-take\" box extends outside image"); }
          }
      }
    else
      { /* Default {TBox} is whole of {I}: */
        for (ax = 0; ax < 2; ax++)
          { o->TBox[ax] = (irange_t){{ 0, o->ISize.c[ax] }}; }
      }
    
    if (argparser_keyword_present(pp, "-pad"))
      { /* Grab and parse the set of paddable edges: */
        ple_parse_edge_set(pp, o->paddable);
      }
    else
      { /* Assume that all edges could be padded: */
        edge_t e; 
        for (e = 0; e < 4; e++) { o->paddable[e] = TRUE; }
      }
      
    argparser_get_keyword(pp, "-minSize");
    o->minSize = ple_parse_size(pp, ple_MAX_SIZE); 
    
    if (argparser_keyword_present(pp, "-midSize"))
      { o->midSize = ple_parse_size(pp, ple_MAX_SIZE);
        int EX0 = o->minSize.c[0], EY0 = o->minSize.c[1];
        int EX1 = o->midSize.c[0], EY1 = o->midSize.c[1];
        /* Check whether {o.midSize} is in the proper range: */
        if ((EX1 < EX0) || (EX1 > 2*EX0) || (EY1 < EY0) || (EY1 > 2*EY0))
          { argparser_error(pp, "min intermediate dimensions are out of range"); }
        /* Check whether {o.midSize} has the same aspect as {o.minSize}: */
        if (EX0*EY1 != EY0*EX1)
          { argparser_error(pp, "min intermediate dimensions have wrong aspect"); }
      } 
    else
      { /* No intermediate series: */
        o->midSize = o->minSize;
      }
      
    argparser_get_keyword(pp, "-maxSize");
    o->maxSize = ple_parse_size(pp, ple_MAX_SIZE); 
    
    o->rotate = argparser_keyword_present(pp, "-rotate");
    
    if (argparser_keyword_present(pp, "-maxScale"))
      { o->maxScale = argparser_get_next_double(pp, 0.0, 1.0); }
    else 
      { o->maxScale = 1.0; }
      
    if (argparser_keyword_present(pp, "-maxStretch"))
      { o->maxStretch = argparser_get_next_double(pp, 0.0, 1.0); }
    else 
      { o->maxStretch = 0.0; }
      
    if (argparser_keyword_present(pp, "-maxPad"))
      { /* Maximum padding factor: */
        o->maxPad = argparser_get_next_double(pp, 0.0, 1.0); }
    else
      { /* No padding allowed: */
        o->maxPad = 0.0;
      }
      
    if (argparser_keyword_present(pp, "-maxCut"))
      { o->maxCut = argparser_get_next_double(pp, 0.0, 1.0); }
    else
      { o->maxCut = 0.5; }
    
    if (argparser_keyword_present(pp, "-debug"))
      { o->debug = argparser_get_next_int(pp, 0, MAXINT); }
    else
      { o->debug = 0; }
    
    /* Check for extraneous arguments: */
    argparser_finish(pp);
        
    return o;
  }
  
i2_t ple_parse_size(argparser_t *pp, int MAX)
  { i2_t sz;
    sz.c[0] = argparser_get_next_int(pp, 0, MAX);
    sz.c[1] = argparser_get_next_int(pp, 0, MAX);
    return sz; 
  }
  
i2_t ple_parse_pos(argparser_t *pp, int MAX)
  { i2_t pos;
    pos.c[0] = argparser_get_next_int(pp, -MAX, MAX);
    pos.c[1] = argparser_get_next_int(pp, -MAX, MAX);
    return pos; 
  }

void ple_parse_edge_set(argparser_t *pp, bool_t es[])
  {
    char *edges = argparser_get_next(pp);
    char *p = edges;
    edge_t e; 
    for (e = 0; e < 4; e++) { es[e] = FALSE; }
    while((*p) != '\000')
      { if ((*p) == 'l')
          { e = edge_L; }
        else if ((*p) == 'r')
          { e = edge_R; }
        else if ((*p) == 't')
          { e = edge_T; }
        else if ((*p) == 'b')
          { e = edge_B; }
        else 
          { argparser_error(pp, "invalid edge letter"); }
        if (es[e]) 
          { argparser_error(pp, "duplicate edge letter"); }
        es[e] = TRUE;
        p++;
      }
  }
