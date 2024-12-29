#define PROG_NAME "pnmhog"
#define PROG_DESC "applies the T-HOG text detector to an image"
#define PROG_VERS "1.0"
/* Last edited on 2024-11-23 06:00:34 by stolfi */

/* Copyright © 2002 by the State University of Campinas (UNICAMP).
** See the copyright, authorship, and warranty notice at end of file.
*/

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  imgc_input_output_coords_HELP " \\\n" \
  "    [ -matrix {D} {TX} {TY}  {DX} {UX} {VX}  {DY} {UY} {VY} \\\n" \
  "    | \\\n" \
  "      -points \\\n" \
  "        {X_IN[1]} {Y_IN[1]} \\\n" \
  "        {X_IN[2]} {Y_IN[2]} \\\n" \
  "        {X_IN[3]} {Y_IN[3]} \\\n" \
  "        {X_IN[4]} {Y_IN[4]} \\\n" \
  "        \\\n" \
  "        {X_OUT[1]} {Y_OUT[1]} \\\n" \
  "        {X_OUT[2]} {Y_OUT[2]} \\\n" \
  "        {X_OUT[3]} {Y_OUT[3]} \\\n" \
  "        {X_OUT[4]} {Y_OUT[4]} \\\n" \
  "    ] \\\n" \
  "    [ -interpolate {INT_ORDER} ] \\\n" \
  "    [ -undef {DEFVAL} | -extend ] \\\n" \
  "    [ -maxval {MV_OUT} ] \\\n" \
  "    [ -isMask {ISMASK} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    [ -debug {XD_OUT} {YD_OUT} ] \\\n" \
  "    [ {PNMFILE_IN} ]"

#define PROG_INFO_DESC \
  "  The program reads the PNM image {PNMFILE_IN} and applies to" \
  " it a specified  projective transformation {M}; so that the" \
  " pixel value at a point {p} of the input image is copied to" \
  " point {M(p)} the output image.\n" \
  "\n" \
  "  The projective map {M} may be specified as a {3×3}" \
  " homogeneous coefficient array (with" \
  " the \"-matrix\" command line argument).  Alternatively, one" \
  " may give four points in the input image and" \
  " four corresponding points in the output image (with \"-points\").  In" \
  " the second case, the matrix {M} is computed using {hr2_pmap_from_four_points}" \
  " in {hr2.h}.\n" \
  "\n" \
  "  In either case, let's represent a generic point {p} of the plane as" \
  " a three-element homogeneous coordinate row vector [{w},{x},{y}], such that" \
  " ({x/w},{y/w}) are the Cartesian coordinats of {p}.  With these conventions," \
  " the vector-matrix product {p·M} (in that order) will" \
  " give the three homogeneous coordinates [{w'},{x'},{y'}]" \
  " of {M(p)}.  If the computed weight {w'} is zero or negative, the point" \
  " {M(p)} is invalid, and the image value at {p} is discarded.\n" \
  "\n" \
  imgc_input_output_coords_intro_INFO("the input image","the output image") "\n" \
  "\n" \
  "  If the argument {PNMFILE_IN} is omitted or is \"-\", the" \
  " program reads the input image from {stdin}.  The output" \
  " image is always written to {stdout}."

#define PROG_INFO_OPTS \
  imgc_parse_input_output_coords_INFO_OPTS( \
    imgc_parse_x_axis_INFO_OPTS_default_pbm, \
    imgc_parse_y_axis_INFO_OPTS_default_pbm, \
    "the input image",
    imgc_parse_center_org_INFO_OPTS_default_zero("-iOrg"), \
    "the output image", \
    imgc_parse_center_org_INFO_OPTS_default_zero("-oOrg"), \
    imgc_parse_size_INFO_OPTS_default_input("-oSize","the input image","the output image") \
  ) "\n" \
  "\n" \
  "  -matrix {D} {TX} {TY}  {DX} {UX} {VX}  {DY} {UY} {VY}\n" \
  "    Specifies the elements of the 3×3 projetive transformation" \
  " matrix, in row by row order.  This argument is mutually" \
  " exclusive with \"-points\".  If \"-matrix\"" \
  " and \"-points\" are both omitted,  the" \
  " projective transformation is trivial (the identity map).\n" \
  "\n" \
  "  -points \\\n" \
  "      {X_IN[1]} {Y_IN[1]} \\\n" \
  "      {X_IN[2]} {Y_IN[2]} \\\n" \
  "      {X_IN[3]} {Y_IN[3]} \\\n" \
  "      {X_IN[4]} {Y_IN[4]}\n" \
  "      \\\n" \
  "      {X_OUT[1]} {Y_OUT[1]} \\\n" \
  "      {X_OUT[2]} {Y_OUT[2]} \\\n" \
  "      {X_OUT[3]} {Y_OUT[3]} \\\n" \
  "      {X_OUT[4]} {Y_OUT[4]}\n" \
  "    This argument states that the projective transformation must" \
  " take each point {(X_IN[i],Y_IN[i])} of the input image to" \
  " the corresponding point {(X_OUT[i],Y_OUT[i])} of the output" \
  " image.  No three of the input-side points may be collinear," \
  " and ditto for the output-side points.  The input and output" \
  " coordinates are relative to the input and output coordinate" \
  " systems, respectively.   This argument is mutually" \
  " exclusive with \"-matrix\".   If \"-points\" and" \
  " \"-matrix\" are both omitted, the" \
  " projective transformation is trivial (the identity map).\n" \
  "\n" \
  "  -interpolate {INT_ORDER} \n" \
  "    This optional argument specifies the way inpt pixels are" \
  " interpolated.  " float_image_transform_interpolation_HELP_INFO "  The default" \
  " is \"-interpolate 0\" (C0 bilinear interpolation).\n" \
  "\n" \
  "  -undef {DEFVAL}\n" \
  "  -extend \n" \
  "    These mutually exclusive optional flags specify the handling" \
  " of source pixels that fall outside the input image's" \
  " domain.  If \"-undef\" is used, any such pixel is assumed to" \
  " have value {DEFVAL}, in a scale from 0 to 1.  If \"-extend\" is used, the input image" \
  " will be implicitly extended to an infinite image, before" \
  " being transformed, by replicating the pixels" \
  " along its borders.  If neither option is specified," \
  " the program assumes \"-undef 0.5\".\n" \
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
  " image statistics.\n" \
  "\n" \
  "  -debug {XD_OUT} {YD_OUT}\n" \
  "    If this option is present, the program prints out" \
  " debugging information about the computation of the" \
  " pixel in column {Xdo} (counting from 0, left to right) and" \
  " row {YD_OUT} (counting from 0, top to bottom) of the output image."

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
  "BUGS\n" \
  "  Sampling is not very scientific; it may blur more" \
  " than necessary, and may not work properly if the map is too warped.\n" \
  "\n" \
  "  The code of ths program is very similar to that of" \
  " \"pnmradist\" and other programs.  That code" \
  " should be shared somehow.\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmscale(1), pnmradist(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created aug/2002 by Jorge Stolfi, IC-UNICAMP as \"pnmgtran\".\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  All changes by J. Stolfi, IC-UNICAMP unless otherwise noted.\n" \
  "\n" \
  "  nov/2006 Rewritten to use sane PNM, argparser, etc.\n" \
  "  jun/2007 Added \"-points\" option.\n" \
  "  jul/2007 Added \"-xAxis\", \"-yAxis\" options.\n" \
  "  aug/2007 Rearranged the points in \"-points\".\n" \
  "  jan/2008 Moved radial distortion to \"pnmradist\".\n" \
  "  jan/2008 Renamed from \"pnmgtran\" to \"pnmprojmap\".\n" \
  "  jan/2008 Added the \"-extend\" option.\n" \
  "  ago/2010 Added the \"-interpolate\" option.\n" \
  "  ago/2010 Added the \"-isMask\" option.\n" \
  "  ago/2023 Updated for changes in {image_input_output_coords.h}.\n" \
  "  ago/2023 Converted {int} to {int32_t}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2002 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

/* We need to set these in order to get {isnan}. What a crock... */
#undef __STRICT_ANSI__
#define _ISOC99_SOURCE 1
#include <math.h>

#include <stdlib.h>
#include <stdio.h>

#include <i2.h>
#include <r2.h>
#include <r2_extra.h>
#include <r2x2.h>
#include <r3x3.h>
#include <hr2.h>

#include <jspnm.h>
#include <interval.h>
#include <ix_reduce.h>
#include <jsfile.h>
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>
#include <image_coords.h>
#include <sample_conv.h>
#include <float_image.h>
#include <float_image_transform.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <argparser.h>

#define MAX_SIZE (32*1024)
  /* A limit on image size, to avoid humongous mallocs. */

typedef struct options_t 
  { char *fname;         /* Input file name. */
    /* User coordinate system axis directions: */
    bool_t yUp;          /* TRUE if the vertical axis points up, FALSE otherwise. */
    bool_t xLeft;        /* TRUE if the horizontal axis points left, FALSE otherwise. */
    /* User coordinates for input image: */
    double iUnit;        /* Unit, in pixels. */
    bool_t iCenter;      /* If TRUE, input origin is center; if FALSE, use {iOrg}. */
    r2_t iOrg;           /* Input origin relative to default origin, if {!iCenter}. */
    /* User coordinates for output image: */
    double oUnit;        /* Unit, in pixels. */
    bool_t oCenter;      /* If TRUE, output origin is center; if FALSE, use {oOrg}. */ 
    r2_t oOrg;           /* Output origin relative to default origin, if {!oCenter}. */
    /* Output image size: */
    double oCols;        /* X size of output image, in user units, or -1 if not given. */
    double oRows;        /* Y size of output image, in user units, or -1 if not given. */
    /* Geometry transformation options: */
    hr2_pmap_t M;        /* The geometric transformation. */
    /* Input image encoding and interpolation options: */
    bool_t isMask;       /* TRUE to interpret samples as in masks. */
    bool_t extend;       /* TRUE extends the the image by row/col replication. */
    float undef;         /* Input image epadding value, if {extend} is false. */
    int32_t interpolate;     /* Interpolation order. */
    /* Output image encoding: */
    uint16_t maxval;     /* Output maxval requested by user, or 0 if not given. */
    /* Debugging options: */
    bool_t verbose;      /* TRUE to print global statistics. */
    i2_t debug;          /* PBM indices of output pixel to debug; or (-1,-1). */
  } options_t;
  /* The transformation conceptually consists of moving the color from
    each point {ip} of the input image's domain (in the user's input
    coordinate system) to the point {op = M.dir(ip)} (in the user's
    output coordinate system). (In practice, we scan each output pixel
    {op} and set it to the input image's value at point {ip =
    M.inv(op)}.  */

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);
options_t *get_options(int32_t argc, char **argv);

void change_coord_systems
  ( hr2_pmap_t *P,    /* (IN/OUT) The projective transformation map. */
    hr2_pmap_t *isys, /* The map from native input to user input coordinates. */
    hr2_pmap_t *osys  /* The map from native output to user output coordinates. */
  );
  /* Assumes that {P} is a projective map from the input domain to the
    output domain, described in terms of the user-defined input and
    output coordinate systems. Modifies the map {P} so that it takes
    and returns points expressed in the native {float_image}
    coordinate systems, according to the given coordinate system options. */

float_image_t *read_image
  ( FILE *rd, 
    int32_t *colsP, 
    int32_t *rowsP, 
    int32_t *chnsP, 
    uint16_t *maxvalP,
    bool_t isMask,
    bool_t verbose
  );
  /* Reads a PBM/PGM/PPM image file from {rd}, converts it to a float
    image with samples in the range [0_1]. Returns the relevant image
    data. If {verbose} is true, prints image statistics to
    {stderr}. */

void write_image
  ( FILE *wr, 
    float_image_t *fim, 
    uint16_t maxval,
    bool_t isMask,
    bool_t verbose
  );
  /* Writes the float image {fim} to {wr} as a PBM/PGM/PPM image file.
    Samples are converted fron the range [0_1]. Returns the relevant
    image data. If {verbose} is true, prints image statistics to
    {stderr}. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    /* Parse command line options: */
    options_t *o = get_options(argc, argv);

    /* Read input image, get dimensions: */
    int32_t chns, iCols, iRows;
    uint16_t maxval_in;
    FILE *rd = open_read(o->fname, o->verbose);
    float_image_t *im_in = read_image(rd, &iCols, &iRows, &chns, &maxval_in, o->isMask, o->verbose);
    
    /* Compute output image size in pixels: */
    int32_t oCols, oRows;
    imgc_compute_output_size_in_pixels(iCols, iRows, o->iUnit, o->oCols, o->oRows, o->oUnit, &oCols, &oRows, MAX_SIZE);
    if (o->verbose) { fprintf(stderr, "output image size (pixels) = %d %d\n", oCols, oRows); }
    
    /* Compute pixel to user coord system map {isys} for input image: */
    hr2_pmap_t isys = imgc_coord_sys_map(o->xLeft, o->yUp, o->iUnit, o->iCenter, &(o->iOrg), iCols, iRows);
    if (o->verbose) { imgc_print_pmap(stderr, "input pixel", "input user", "isys", &(isys)); }

    /* Compute pixel to user coord system map {osys} for the output image: */
    hr2_pmap_t osys = imgc_coord_sys_map(o->xLeft, o->yUp, o->oUnit, o->oCenter, &(o->oOrg), oCols, oRows);
    if (o->verbose) { imgc_print_pmap(stderr, "output pixel", "output user", "osys", &(osys)); }

    /* Build a modified projective map {N} acting on {float_image} coords: */
    if (o->verbose) { imgc_print_pmap(stderr, "input user", "output user", "M", &(o->M)); }
    hr2_pmap_t N = o->M;
    change_coord_systems(&N, &isys, &osys);
    if (o->verbose) { imgc_print_pmap(stderr, "input pixel", "output pixel", "N", &N); }
   
    /* Check whether projective map is the identity: */
    bool_t ident_proj = hr2_pmap_is_identity(&N);
    
    /* Allocate output image: */
    float_image_t *im_ot = float_image_new(chns, oCols, oRows);
    
    /* Map center of debugging pixel from PBM to native output coords {dbp}: */
    r2_t dbp;
    bool_t dbp_defined;
    if ((o->debug.c[0] < 0) || (o->debug.c[1] < 0))
      { dbp = (r2_t){{ NAN, NAN }};
        dbp_defined = FALSE;
      }
    else
      { dbp.c[0] = o->debug.c[0];
        dbp.c[1] = oRows - o->debug.c[1];
        r2_t odbp;
        r2_map_projective(&dbp, &(osys.dir), &odbp, NULL);
        dbp_defined = TRUE;
        if (o->verbose) 
          { fprintf(stderr, "watching output point with coordinates:\n");
            fprintf(stderr, " (%7d,%7d) PBM (Y down)\n", o->debug.c[0], o->debug.c[1]);
            fprintf(stderr, " (%7.1f,%7.1f) libimg (Y up)\n", dbp.c[0], dbp.c[1]);
            fprintf(stderr, " (%7.1f,%7.1f) user\n", odbp.c[0], odbp.c[1]);
          }
      }
    
    /* Compute output image from input image: */
    
    auto void map_point(r2_t *pP, r2x2_t *JP);
      /* Maps a point {*pP} of the OUTPUT image (in the {float_image}
        native coordinates, with (0,0) at bottom left) to the
        corresponding point of the INPUT image, and stored that into {*pP}.  
        Also multplies {*JP} by the jacobian. If the point
        falls outside the input image's domain, and {o->extend} is
        false, sets {*pP} to {(NAN,NAN)}. */
        
    auto bool_t debug_point(r2_t *pP);
      /* True if {*pP} is the watched pixel. */
      
    ix_reduce_mode_t red = ix_reduce_mode_SINGLE;
    float_image_transform_all(im_in, red, &map_point, o->undef, TRUE, o->interpolate, &debug_point, im_ot);
    
    /* Choose output maxval: */
    uint16_t maxval_ot = (o->maxval > 0 ? o->maxval : maxval_in);
    if (maxval_ot < 255) { maxval_ot = 255; }

    write_image(stdout, im_ot, maxval_ot, o->isMask, o->verbose);

    if (o->verbose) { fprintf(stderr, "done.\n"); }
    return 0;
    
    bool_t debug_point(r2_t *pP)
      { return dbp_defined &&
          (fabs(pP->c[0] - dbp.c[0]) <= 0.501) && 
          (fabs(pP->c[1] - dbp.c[1]) <= 0.501);
      }
    
    void map_point(r2_t *pP, r2x2_t *JP)
      {
        /* Apply the inverse projective map: */
        if (! ident_proj) { r2_map_projective(pP, &(N.inv), pP, JP); }
          
        if (! o->extend)
          { /* Check domain: */
            bool_t invalid = FALSE;
            invalid |= ((pP->c[0] < 0) || (pP->c[0] > iCols));
            invalid |= ((pP->c[1] < 0) || (pP->c[1] > iRows));
            if (invalid) { pP->c[0] = pP->c[1] = NAN; }
          }
      }
  }

void change_coord_systems
  ( hr2_pmap_t *P,
    hr2_pmap_t *isys, /* The map from native input to user input coordinates. */
    hr2_pmap_t *osys  /* The map from native output to user output coordinates. */
  )
  {
    /* Replace {P} by {(isys)*P*(osys^-1): */
    (*P) = hr2_pmap_compose(isys, P);
    hr2_pmap_t syso = hr2_pmap_inv(osys);
    (*P) = hr2_pmap_compose(P, &syso);
  }

#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

options_t *get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
     
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    /* Set defaults for input and output coord systems: */
    o->xLeft = FALSE;
    o->yUp = FALSE;
    o->iCenter = FALSE;
    o->iOrg = (r2_t){{ 0.0, 0.0 }};
    o->oCenter = FALSE;
    o->oOrg = (r2_t){{ 0.0, 0.0 }};
    
    /* The default output size depends on the input image size, so leave {-1}: */
    o->oCols = -1.0; 
    o->oRows = -1.0;
    
    /* Parse input and output coord systems, and output size: */
    imgc_parse_input_output_coords_args
      ( pp, &(o->xLeft), &(o->yUp), 
        &(o->iUnit), &(o->iCenter), &(o->iOrg), 
        &(o->oUnit), &(o->oCenter), &(o->oOrg), 
        &(o->oCols), &(o->oRows)
      );
    
    o->M = argparser_get_proj_map(pp);
    
    if (argparser_keyword_present(pp, "-interpolate"))
      { o->interpolate = (int32_t)argparser_get_next_int(pp, 0, 1); }
    else
      { /* The default depends on the input image: */
        o->interpolate = 0;
      }
    
    if (argparser_keyword_present(pp, "-extend"))
      { o->extend = TRUE; o->undef = 0.5; }
    else if (argparser_keyword_present(pp, "-undef"))
      { o->extend = FALSE;
        o->undef = (float)argparser_get_next_double(pp, 0.0, 1.0);
      }
    else
      { o->extend = FALSE; o->undef = 0.5; }
    
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (uint16_t)argparser_get_next_int(pp, 1, PNM_MAX_SAMPLE); }
    else
      { /* The default depends on the input image: */
        o->maxval = 0;
      }

    if (argparser_keyword_present(pp, "-isMask"))
      { o->isMask = argparser_get_next_bool(pp); }
    else
      { o->isMask = FALSE; }

    o->verbose = argparser_keyword_present(pp, "-verbose");

    if (argparser_keyword_present(pp, "-debug"))
      { o->debug.c[0] = (int32_t)argparser_get_next_int(pp, -1, MAX_SIZE);
        o->debug.c[1] = (int32_t)argparser_get_next_int(pp, -1, MAX_SIZE);
      }
    else 
      { /* The defaults depend on the input image: */
        o->debug = (i2_t){{ -1, -1 }};
      }
    
    /* Parse optional input file name: */
    argparser_skip_parsed(pp);
    if (argparser_next(pp) != NULL) 
      { o->fname = argparser_get_next(pp); }
    else
      { o->fname = "-"; }

    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }

float_image_t *read_image
  ( FILE *rd, 
    int32_t *colsP, 
    int32_t *rowsP, 
    int32_t *chnsP, 
    uint16_t *maxvalP,
    bool_t isMask,
    bool_t verbose
  )
  { uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, TRUE, verbose);
    (*colsP) = pim->cols;
    (*rowsP) = pim->rows;
    (*chnsP) = pim->chns;
    (*maxvalP) = pim->maxval;
    uint16_image_free(pim);
    return fim;
  }

void write_image
  ( FILE *wr, 
    float_image_t *fim, 
    uint16_t maxval,
    bool_t isMask,
    bool_t verbose
  )
  { int32_t chns = (int32_t)fim->sz[0];
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, chns, NULL, NULL, NULL, maxval, TRUE, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    uint16_image_free(pim);
  }
 
