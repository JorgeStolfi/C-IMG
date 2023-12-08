/* Last edited on 2023-10-10 19:10:44 by stolfi */

#define PROG_NAME "pnmpralign"
#define PROG_DESC "Finds a projective map that aligns two images"
#define PROG_VERS "1.0"

#define pnmpralign_C_COPYRIGHT \
  "Copyright © 2017 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    {IMG1} \\\n" \
  "    {IMG2} \\\n" \
  "    " imgc_parse_x_axis_HELP " \\\n" \
  "    " imgc_parse_y_axis_HELP " \\\n" \
  "    " imgc_parse_center_org_HELP("-center","-org","") " \\\n" \
  "    { -feature {X1k} {Y1k} {R1k} {X2k} {Y2k} {R2k} {Wk}  }.. \\\n" \
  "    [ -undef {DEFVAL} | -extend ] \\\n" \
  "    [ -isMask {ISMASK} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO_DESC \
  "  The program reads two image files {IMG1} and {IMG2} and finds a projective map {PM}" \
  " that transforms one to the other.  Both files should be in a Netpbm file" \
  " format (\".pbm\", \".pgm\", or \".ppm\").\n" 
  "\n" \
  "  More precisely, looks for a projective map {QM} of the plane that" 
  " maps {IMG1} to a /reference plane/ an then the reference plane" 
  " to {IMG2}.  The map {QM} is chosen so that as to minimize the" 
  " total squared mismatch {U2(QM)} between {QM(IMG1)}" 
  " and {QM^{-1}(IMG2)}.  The metric {U2(QM)} is defined as" 
  " the sum of the weighted squared difference between {IMG1(x1,y1)}" 
  " and {IMG2(x2,y2)}, where {(x1,y1)=QM^{-1}(x,y)}" 
  " and {(x2,y2)=QM(x,y)}, over a certain set {SMP} of" 
  " sampling points {(x,y)} on the reference plane.  Then" 
  " it returns the map {PM=QM^2}.\n" 
  "\n" \
  "  The image values are interpolated when computing the" 
  " mismatch.  The sampling points and their weights" 
  " are specified by the \"-grid\" or \"-feature\" command line arguments.\n" 
  "\n" \
  "  The program ignores the extensions of the input filenames, and" \
  " interprets the contents according to the 16-bit \"magic number\" at" \
  " the start of the file.  All sample values are converted from integers" \
  " in {0..MAXVAL} to floating-point values" \
  " between 0 and 1; where {MAXVAL} is the maximum sample value as" \
  " declared in the input file.   Note that in PBM (binary) files" \
  " the encoding is reversed (0-->1, 1-->0).\n\n" \
  "\n" \
  "  " argparser_proj_map_INFO ".\n" \
  "\n" \
  "  "  "\n" \
  "\n" \
  "  " imgc_user_axes_INFO "" \
  "  " imgc_pixel_axes_INFO "" \
  "  " imgc_pixel_centers_INFO "" \
  "  " imgc_origin_INFO("-center","-org","both images") ""

#define PROG_INFO_OPTS \
  imgc_parse_x_axis_INFO_OPTS(imgc_parse_x_axis_INFO_OPTS_default_pbm) "" \
  "  This parameter affects the" \
  " interpretation of all X coordinates in the arguments, including" \
  " {CX1}, {CX2}, {X1[k]}, {X1[k]}, and the matrix" \
  " coefficients.\n" \
  "\n" \
  imgc_parse_y_axis_INFO_OPTS(imgc_parse_y_axis_INFO_OPTS_default_pbm) "" \
  "  This parameter affects the interpretation of all Y" \
  " coordinates in the arguments, including" \
  " {CY1}, {CY2}, {Y1[k]}, {Y2[k]}, and the matrix" \
  " coefficients.\n" \
  "\n" \
  imgc_parse_center_org_INFO_OPTS("-center","-org","",imgc_parse_center_org_INFO_OPTS_default_zero("-org")) "\n" \
  "\n" \
  "  -grid {NX} {NY} \n" \
  "    This option specifies that the  sampling grid" 
  " has {NX} columns and {NY} rows of sampling points, placed" 
  " so that {(x1,y1)} and {(x2,y2)} cover most of {IMG1}" 
  " and {IMG2}, respectively.  The weight {W(x,y)} used for" 
  " each sampling point {(x,y)} decays smoothly" 
  " from 1 to 0 as {(x,y)} gets close to the edges of the grid.\n" \
  "\n" \
  "  -feature {X1k} {Y1k} {R1k} {X2k} {Y2k} {R2k} {Wk} \n" \
  "    This option specifies that one term of the mismatch function" 
  " will the point difference between" 
  " so that {(x1,y1)} and {(x2,y2)} cover most of {IMG1}" 
  " and {IMG2}, respectively.  The weight {W(x,y)} used for" 
  " each sampling point {(x,y)} decays smoothly" 
  " from 1 to 0 as {(x,y)} gets close to the edges of the grid.\n" \
  "\n" \
  argparser_proj_map_INFO_OPTS "\n" \
  "\n" \
  "  -interpolate {INT_ORDER} \n" \
  "    This optional argument specifies the way inpt pixels are" \
  " interpolated.  " float_image_transform_interpolation_INFO_OPTS "  The default" \
  " is \"-interpolate 0\" (C0 bilinear interpolation).\n" \
  "\n" \
  "  -undef {DEFVAL}\n" \
  "  -extend \n" \
  "    These mutually exclusive optional flags specify the handling" \
  " of source pixels that fall outside the input image's" \
  " domain.  If \"-undef\" is used, any such pixel is assumed to" \
  " have value {DEFVAL}, in a scale" 
  " from 0 to 1.  If \"-extend\" is used, the input image" \
  " will be implicitly extended to an infinite image, before" \
  " being transformed, by replicating the pixels" \
  " along its borders.  If neither option is specified," \
  " the program assumes \"-undef 0.5\".\n" \
  "\n" \
  imgc_parse_output_size_INFO_OPTS "  If omitted, the output image" \
  " will have the same size as the input one.\n" \
  "\n" \
  "  -maxval {MV_OUT}\n" \
  "    Specifies {MV_OUT} as the maximum sample value for the" \
  " output image.  It must be an integer between 255 and 65535," \
  " inclusive. If not specified, it is set to 255 or to the" \
  " input image's {maxval}, whichever is larger.\n" \
  "\n" \
  "  -isMask {ISMASK}\n" \
  "    This optional Boolean argument modifies the interpretation" \
  " of integer sample values, especially 0 and {MAXVAL}, in the" \
  " input and output files (see {sample_conv.h}).  If {ISMASK} is true (\"T\" or 1)," \
  " " sample_conv_0_1_isMask_true_INFO "  If {ISMASK} is false (\"F\" or 0)," \
  " " sample_conv_0_1_isMask_false_INFO "  The default is \"-isMask F\".\n" \
  "\n" \
  "  -verbose\n" \
  "    If this option is present, the program prints out" \
  " global debugging information, such as input and output" \
  " image statistics."

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
  "INPUT FILES\n" \
  "  " float_image_read_gen_INFO "\n" \
  
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_INFO_OPTS "\n" \
  "\n" \
  "BUGS\n" \
  "  The code of ths program is very similar to that of" \
  " \"pnmprojmap\" and other programs.  That code" \
  " should be shared somehow.\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmprojmap(1), pnmradist(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2017-03-03 by Jorge Stolfi, IC-UNICAMP from bits of \"pnmprojmap\", \"fni_view\", etc.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  mar/2017 Created.  J. Stolfi, IC-UNICAMP.\n" \
  "  jun/2017 Worked some more.  J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmpralign_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <r3.h>
#include <r3x3.h>
#include <hr3.h>
#include <float_image.h>
#include <float_image_read_gen.h>
#include <image_coords.h>

/* COMMAND-LINE OPTIONS */

typedef struct coord_sys_t 
  { /* Global coordinate system and encoding options: */
    bool_t yUp;          /* TRUE if the vertical axis points up, FALSE otherwise. */
    bool_t xLeft;        /* TRUE if the horizontal axis points left, FALSE otherwise. */
    /* Image-specific coordinate system options: */
    bool_t iCenter;      /* If TRUE, input origin is center; if FALSE, use {iOrg}. */
    r2_t iOrg;           /* Input origin relative to default origin, if {!iCenter}. */
    bool_t oCenter;      /* If TRUE, output origin is center; if FALSE, use {oOrg}. */ 
    r2_t oOrg;           /* Output origin relative to default origin, if {!oCenter}. */
  } coord_options_t;
  /* Image coordinate system specs. */

typedef struct options_t
  { char *fnameA;        /* Image A file name. */
    char *fnameB;        /* Image B file name. */
    /* Input image coordinates, encoding, and interpolation options: */
    coord_sys_t sys;     /* Image coordinate systems. */
    bool_t isMask;       /* TRUE to interpret samples as in masks. */
    bool_t extend;       /* TRUE extends the the image by row/col replication. */
    float undef;         /* Input image epadding value, if {extend} is false. */
    int interpolate;     /* Interpolation order. */
    /* Debugging options: */
    bool_t verbose;      /* TRUE to print global statistics. */
  } options_t;

/* INTERNAL PROTOTYPES */

options_t *pnmpralign_parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int main(int argc,char** argv);

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  {
    /* Parse command line: */
    options_t *o = pnmpralign_parse_options(argc, argv);
    
    
    
    return 0;
}

options_t *pnmpralign_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    /* Parse keyword parameters: */
    
    /* Parse the boolean option {op1}: */
    o->op1 = argparser_keyword_present(pp, "-op1");
    
    /* Parse the string option {op2}: */
    if (argparser_keyword_present(pp, "-op2"))
      { o->op2 = argparser_get_next(pp); }
    else
      { o->op2 = "NONE"; }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    if (argparser_next(pp) != NULL)
      { o->infile = argparser_get_next(pp); }
    else
      { o->infile = "-"; }

    if (argparser_next(pp) != NULL)
      { o->outfile = argparser_get_next(pp); }
    else
      { o->outfile = "-"; }

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

