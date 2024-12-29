#define PROG_NAME "fni_crop"
#define PROG_DESC "crops a FNI (float image) file"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-25 07:15:09 by stolfi */

#define PROG_C_COPYRIGHT "Copyright © 2024 Jorge Stolfi, UNICAMP. Run \"" PROG_NAME " -info\" for details"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <argparser.h>
#include <jsfile.h>
#include <float_image.h>
#include <image_coords.h>
#include <limits.h>
#include <math.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    [ -top {NUM} ] \\\n" \
  "    [ -left {NUM} ] \\\n" \
  "    [ -width {NUM} ] \\\n" \
  "    [ -height {NUM} ] \\\n" \
  "    [ -right {NUM} ] \\\n" \
  "    [ -bottom {NUM} ] \\\n" \
  "    [ -bottom {NUM} ] \\\n" \
  "    " imgc_parse_y_axis_HELP " \\\n" \
  "    [ -verbose ] \\\n" \
  "    [ -in {IN_FILE} ] \\\n" \
  "    [ -out {OUT_FILE} ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION" \
  "  This program reads the file {IN_FILE} and writes a rectangular" \
  " sub-image of it to {OUTPUT_FILE}.  Both files are assumed to be" \
  " in FNI (float image) format.  The sub-image may extend outside" \
  " the original image, in which case the non-existing pixels are set to {NAN}.\n" \
  "\n" \
  "  The \"-yAxis\" option determines whether row 0 of the" \
  " FNI image is the top or bottom row.\n" \
  "\n" \
  "OPTIONS" \
  "  -width {WIDTH}\n" \
  "  -height {HEIGHT}\n" \
  "    These optional arguments specify the width and height of the sub-image." \
  "\n" \
  "  -top {TOP}\n" \
  "  -left {LEFT}\n" \
  "  -right {RIGHT}\n" \
  "  -bottom {BOTTOM}\n" \
  "    These optional arguments ask for the removal of the specidied rows or columns" \
  " at the specified side of the image.  The arguments may be negative, in which" \
  " case rows or columns of {NAN} pixels are added instead of trimmed.\n" \
  "\n" \
  "    The sum of {LEFT}, {RIGHT}, and {WIDTH} must be equal to the input image" \
  " width.  If only one of them" \
  " is omitted, it is assumed to be whatever is" \
  " needed to make up to the input image width.   If two or more are omitted, then" \
  " missing {LEFT} and/or {RIGHT} (in that order) are assumed to be zero.  In any" \
  " case, the final {WIDTH} must be non-negative.\n" \
  "\n" \
  "    Similar rules apply to {TOP}, {BOTTOM}, and {HEIGHT}.\n" \
  "\n" \
  imgc_parse_y_axis_INFO_OPTS(imgc_parse_y_axis_INFO_OPTS_default_pbm) "\n" \
  "\n" \
  "  -in {IN_FILE}\n" \
  "    This optional argument specifies the name of the input file.  If" \
  " omitted, the program reads from {stdin}." \
  "\n" \
  "  -verbose\n" \
  "    If this optional argument is present, the program prints some" \
  " informative messages to {stderr}." \
  "\n" \
  "  -out {OUT_FILE}\n" \
  "    This optional argument specifies the name of the output file.  If" \
  " omitted, the program writes to {stdout}." \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  Santorini.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created oct/2024 by Jorge Stolfi (UNICAMP).\n" \
  "  Based on {fni_cut} created in 2010 or earlier by Rafael Saracchini.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2024-12-24 J.Stolfi Renamed; added width/height; negative values.\n" \
  "  2024-12-24 J.Stolfi Added \"-flipY\", \"-verbose\".\n" \
  "  2024-12-24 J.Stolfi Added \"-yAxis\".\n"  \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

typedef struct options_t
  { int32_t width;
    int32_t height;
    int32_t top;
    int32_t left;
    int32_t bottom;
    int32_t right;
    bool_t yaxis_up;
    bool_t verbose;
    char* in;
    char* out;
    bool_t yUp;
  } options_t;

void fix_defaults(int32_t N, int32_t *lo_P, int32_t *hi_P, int32_t *sz_P, bool_t verbose);
  /* Takes an image dimension {N} (count of either columns or rows) and 
    fixes the low-end trim count {*lo_P}, the high-end trim count {*hi_P},
    and the remaining size {*sz_P}, providing suitable defaults for unspecified
    ones.  Assumes that a parameter is unspecified if and only if
    ]its value is {INT32_MAX}. */

float_image_t* read_image (char* filename, bool_t verbose);
  /* Read the file {filename}, or {stdin} if  {filename} is {NULL}. */

float_image_t* cut_image (float_image_t* img_in, int32_t left, int32_t right, bool_t yUp, int32_t top, int32_t bottom);
  /* Applies the cropping specified by {left,right,top,bottom} to {img_in}.
    If {yUp} is true, assumes that the bottom of the image is row 0, otherwise
    assumed that it is row {NY-1}. */

void write_image (char* filename, float_image_t* img, bool_t verbose);
  /* Writes and FNI file, or to {stdout} if filename is {NULL}. */

options_t* parse_args (int32_t argc, char** argv);
  /* Parses the command line arguments. */

 int32_t main(int32_t argc, char** argv); 
  
/* IMPLEMENTATONS */

int32_t main(int32_t argc, char** argv)
  {
    options_t* o = parse_args(argc,argv);
    float_image_t* img_in = read_image(o->in, o->verbose);

    int32_t NC_in, NX_in, NY_in;
    float_image_get_size(img_in, &NC_in, &NX_in, &NY_in);

    fix_defaults(NX_in, &(o->left), &(o->right), &(o->width), o->verbose);
    if (o->verbose) 
      { fprintf(stderr, "input width = %d - trimming/adding columns", NX_in);
        fprintf(stderr, " %d at left, %d at right, output width %d\n", o->left, o->right, o->width);
        int32_t ini_in = o->left;
        int32_t fin_in = NX_in - 1 - o->right;
        fprintf(stderr, " output cols {0..%d} are input cols {%d..%d}\n", o->width-1, ini_in, fin_in);
      }

    fix_defaults(NY_in, &(o->top), &(o->bottom), &(o->height), o->verbose);
    if (o->verbose) 
      { fprintf(stderr, "input height = %d - trimming/adding rows", NY_in);
        fprintf(stderr, " %d at top, %d at bottom, output height %d\n", o->top, o->bottom, o->height);
        int32_t ini_in = (o->yaxis_up ? o->bottom : o->top);
        int32_t fin_in = NY_in - 1 - (o->yaxis_up ? o->top : o->bottom);
        fprintf(stderr, " output rows {0..%d} are input rows {%d..%d}\n", o->width-1, ini_in, fin_in);
      }

    float_image_t* img_ot = cut_image(img_in, o->left, o->right, o->yaxis_up, o->top, o->bottom);
    write_image(o->out,img_ot, o->verbose);
    float_image_free(img_in);
    float_image_free(img_ot);
    return 0;
  }
  
void fix_defaults(int32_t N, int32_t *lo_P, int32_t *hi_P, int32_t *sz_P, bool_t verbose)
  { int32_t lo = (*lo_P);
    int32_t hi = (*hi_P);
    int32_t sz = (*sz_P);
    
    if (sz != INT32_MAX)
      { /* Final size is spefcified. Provide defaults for {lo} and {hi}, in that order: */
        if (lo != INT32_MAX)
          { if (hi == INT32_MAX) { hi = N - lo - sz; } }
        else if (hi != INT32_MAX)
          { lo = N - sz - hi; }
        else
          { lo = 0; hi = N - sz; }
      }
    else
      { /* Final size is unspecified, compute from {lo,hi}: */
        if (lo == INT32_MAX) { lo = 0; }
        if (hi == INT32_MAX) { hi = 0; }
        sz = N - lo - hi;
      }
    /* Now all three must be defined: */
    if (verbose) fprintf(stderr, "  N = %d  lo = %d  sz = %d  hi = %d  tot = %d\n", N, lo, sz, hi, lo+sz+hi);
    demand(N == lo + sz + hi, "trim amounts and output size don't add to input image size");
    demand(sz >= 0, "invalid (negative) output image size");
    (*lo_P) = lo;
    (*hi_P) = hi;
    (*sz_P) = sz;
  }

float_image_t* cut_image (float_image_t* img_in, int32_t left, int32_t right, bool_t yUp, int32_t top, int32_t bottom) 
  { int32_t NX_in,NY_in,NC_in;
    float_image_get_size(img_in, &NC_in, &NX_in, &NY_in);
    
    int32_t NC_ot,NX_ot,NY_ot;
    NC_ot = NC_in;
    NX_ot = NX_in - left - right; demand(NX_ot >= 0, "excessive {left,right} trim");
    NY_ot = NY_in - top - bottom; demand(NY_ot >= 0, "excessive {top,bottom} trim");
    float_image_t* img_ot = float_image_new(NC_ot,NX_ot,NY_ot);

    for(int32_t ox = 0; ox < NX_ot; ox++)
      { for(int32_t oy = 0; oy < NY_ot; oy++)
          { for(int32_t c = 0; c < NC_ot; c++)
              { int32_t ix = ox + left;
                int32_t iy = oy + (yUp ? bottom : top);
                float val;
                if ((ix < 0) || (ix >= NX_in) || (iy < 0) || (iy >= NY_in))
                  { val = NAN; }
                else
                  { val = float_image_get_sample(img_in, c, ix, iy); }
                float_image_set_sample(img_ot, c, ox, oy, val);
              }
          }
      }
    return img_ot;
  }

options_t* parse_args(int32_t argc, char** argv)
  {
    options_t* o = (options_t*)malloc(sizeof(options_t));
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    o->width = INT32_MAX;
    o->height = INT32_MAX;
    o->left = INT32_MAX;
    o->right = INT32_MAX;
    o->top = INT32_MAX;
    o->bottom = INT32_MAX;

    if (argparser_keyword_present(pp, "-left"))
      { o->left = (int32_t)argparser_get_next_int(pp, INT32_MIN+2, INT32_MAX-1); }
    if (argparser_keyword_present(pp, "-right"))
      { o->right = (int32_t)argparser_get_next_int(pp, INT32_MIN+2, INT32_MAX-1); }
    if (argparser_keyword_present(pp, "-top"))
      { o->top = (int32_t)argparser_get_next_int(pp, INT32_MIN+2, INT32_MAX-1); }
    if (argparser_keyword_present(pp, "-bottom"))
      { o->bottom = (int32_t)argparser_get_next_int(pp, INT32_MIN+2, INT32_MAX-1); }

    if (argparser_keyword_present(pp, "-width"))
      { o->width = (int32_t)argparser_get_next_int(pp, INT32_MIN+2, INT32_MAX-1); }
    if (argparser_keyword_present(pp, "-height"))
      { o->height = (int32_t)argparser_get_next_int(pp, INT32_MIN+2, INT32_MAX-1); }

    o->yUp = FALSE;
    imgc_parse_y_axis(pp, &(o->yUp));
     
    o->verbose =  argparser_keyword_present(pp, "-verbose");

    o->in = NULL;
    if (argparser_keyword_present(pp, "-in"))
      { o->in =  argparser_get_next(pp); }

    o->out = NULL;
    if (argparser_keyword_present(pp, "-out"))
      { o->out =  argparser_get_next(pp); }

    argparser_finish(pp);
    return o;
  }

float_image_t* read_image(char* filename, bool_t verbose)
  {
    FILE *rd;
    if (filename == NULL)
      { if (verbose) { fprintf(stderr, "reading image from {stdin} ...\n"); }
        rd = stdin;
      }
    else
      { rd = open_read(filename, verbose); };
    float_image_t* img_in = float_image_read(rd);
    if(filename != NULL) { fclose(rd); }
    return img_in;
  }

void write_image(char* filename, float_image_t* img, bool_t verbose)
  { FILE* wr;
    if (filename == NULL)
      { if (verbose) { fprintf(stderr, "writing image to {stdout} ...\n"); }
        wr = stdout;
      }
    else
      { wr = open_write(filename, verbose); };
    float_image_write(wr,img);
    if(filename != NULL) { fclose(wr); } else { fflush(wr); }
  }
