#define PROG_NAME "pnmradist"
#define PROG_DESC "applies radial distortion to a pbm/ppm/pgm file"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-23 05:59:25 by stolfi */

/* Copyright © 2002 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  imgc_input_output_coords_HELP " \\\n" \
  "    -kappa {KAPPA} \\\n" \
  "    -pixelSize {HX} {HY} \\\n" \
  "    [ -interpolate {INT_ORDER} ] \\\n" \
  "    [ -extend ] \\\n" \
  "    [ -maxval {MV_OUT} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    [ -debug {XD_OUT} {YD_OUT} ] \\\n" \
  "    [ {PNMFILE_IN} ]"

#define PROG_INFO_DESC \
  "  The program reads the PNM image {PNMFILE_IN} and applies to" \
  " it a radial distortion with parameter {KAPPA}.\n" \
  "\n" \
  "  The radial correction is applied" \
  " relative to the specified origin of the input image.\n" \
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
    "the input image", \
    imgc_parse_center_org_INFO_OPTS_default_center("-iCenter"), \
    "the output image", \
    imgc_parse_center_org_INFO_OPTS_default_center("-oCenter"), \
    imgc_parse_size_INFO_OPTS_default_input("-oSize","the input image","the output image") \
  ) "\n" \
  "\n" \
  "  -kappa {KAPPA}\n" \
  "    This argument specifies the amount of radial" \
  " distortion formula; its units are {mm^{-2}}, measured" \
  " on the sensor plane.  Positive" \
  " values of {KAPPA} give a pincushion" \
  " distortion, negative values give a barrel distortion. If" \
  " {KAPPA} is zero, there is no correction.\n" \
  "\n" \
  "  -pixelSize {HX} {HY}\n" \
  "    This argument defines the nominal width and" \
  " height of each pixel sensor in millimeters.  Note that" \
  " these values affoec the interpretation of the parameter {KAPPA}.\n" \
  "\n" \
  "  -interpolate {INT_ORDER} \n" \
  "    This optional argument specifies the way inpt pixels are" \
  " interpolated.  " float_image_transform_interpolation_HELP_INFO "  The default" \
  " is \"-interpolate 0\" (C0 bilinear interpolation).\n" \
  "\n" \
  "  -extend \n" \
  "    This optional flag specifies that the input image is" \
  " to be implicitly extended to a n infinite image, before" \
  " being subjected to the map, by replicating the pixels" \
  " along its borders.  If this option is omitted, pixels" \
  " outside the input image's domain are assumed to be undefined." \
  "\n" \
  "  -maxval {MV_OUT}\n" \
  "    Specifies {MV_OUT} as the maximum sample value for the" \
  " output image.  It must be an integer between 255 and 65535," \
  " inclusive. If not specified, it is set to 255 or to the" \
  " input image's {maxval}, whichever is larger.\n" \
  "\n" \
  "  -verbose\n" \
  "    If this option is present, the program prints out" \
  " global debugging information, such as input and output" \
  " image statistics.\n" \
  "\n" \
  "  -debug {XD_OUT} {YD_OUT}\n" \
  "    If this option is present, the program prints out" \
  " debugging information about the computation of the" \
  " pixel in column {XD_OUT} (counting from 0, left to right) and" \
  " row {YD_OUT} (counting from 0, top to bottom) of the output" \
  " image, regardless of the coordinate system being used."

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
  " \"pnmprojmap\" and other programs.  That code" \
  " should be chared somehow.\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmscale(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created jan/2009 by Jorge Stolfi, IC-UNICAMP, from {pnmtran.c}.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  All changes by J. Stolfi, IC-UNICAMP unless otherwise noted.\n" \
  "\n" \
  "  jan/2008 Added the \"-extend\" option.\n" \
  "  ago/2010 Added the \"-interpolate\" option.\n" \
  "  ago/2023 Updated for changes in {image_input_output_coords.h}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2009 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

/* We need to set these in order to get {isnan}. What a crock... */
#undef __STRICT_ANSI__
#define _ISOC99_SOURCE 1
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <r2.h>
#include <r2_extra.h>
#include <i2.h>
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
#include <float_image.h>
#include <float_image_transform.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <argparser.h>

#define MAX_SIZE (32*1024)
  /* A limit on image size, to avoid humongous mallocs. */

typedef struct options_t 
  { char *fname;         /* Input file name. */
    /* Global coordinate system options: */
    bool_t yUp;          /* TRUE if the vertical axis points up, FALSE otherwise. */
    bool_t xLeft;        /* TRUE if the horizontal axis points left, FALSE otherwise. */
    /* User coordinates for input image: */
    double iUnit;          /* Unit, in pixels. */
    bool_t iCenter;      /* If TRUE, input origin is center; if FALSE, use {iOrg}. */
    r2_t iOrg;           /* Input origin relative to default origin, if {!iCenter}. */
    /* User coordinates for output image: */
    double oUnit;          /* Unit, in pixels. */
    bool_t oCenter;      /* If TRUE, output origin is center; if FALSE, use {oOrg}. */ 
    r2_t oOrg;           /* Output origin relative to default origin, if {!oCenter}. */
    /* Output image attributes: */
    double oCols;        /* X size of output image, in user units, or -1 if not given. */
    double oRows;        /* Y size of output image, in user units, or -1 if not given. */
    /* Geometry transformation options: */
    double kappa;        /* Radial deformation parameter (0 = none). */
    r2_t pixelSize;      /* Pixel dimensions in {mm}. */
    bool_t extend;       /* TRUE extends the the image by row/col replication. */
    int32_t interpolate;     /* Interpolation order. */
    /* Output image encoding: */
    uint16_t maxval;       /* Output maxval requested by user, or 0 if not given. */
    /* Debugging options: */
    bool_t verbose;      /* TRUE to print global statistics. */
    i2_t debug;          /* PBM indices of output pixel to debug; or (-1,-1). */
  } options_t;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);
options_t *get_options(int32_t argc, char **argv);

float_image_t *read_image
  ( FILE *rd, 
    int32_t *colsP, 
    int32_t *rowsP, 
    int32_t *chnsP, 
    uint16_t *maxvalP,
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
    float_image_t *im_in = read_image(rd, &iCols, &iRows, &chns, &maxval_in, o->verbose);
     
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
        corresponding point of the INPUT image, and stores it into {*pP}. 
        Also multplies {*JP} by the jacobian. If the point
        falls outside the input image's domain, and {o->extend} is
        false, sets {*pP} to {(NAN,NAN)}. */
        
    auto bool_t debug_point(r2_t *pP);
      /* True if {*pP} is the watched pixel. */
      
    if (o->verbose) { fprintf(stderr, "transforming the image ...\n"); }
    ix_reduce_mode_t red = ix_reduce_mode_SINGLE;
    float_image_transform_all(im_in, red, &map_point, 0.5, TRUE, o->interpolate, &debug_point, im_ot);
    
    /* Choose the output maxval: */
    uint16_t maxval_ot = (o->maxval > 0 ? o->maxval : maxval_in);
    if (maxval_ot < 255) { maxval_ot = 255; }

    write_image(stdout, im_ot, maxval_ot, o->verbose);

    if (o->verbose) { fprintf(stderr, "done.\n"); }
    return 0;
    
    bool_t debug_point(r2_t *pP)
      { return dbp_defined &&
          (fabs(pP->c[0] - dbp.c[0]) <= 0.501) && 
          (fabs(pP->c[1] - dbp.c[1]) <= 0.501);
      }
    
    void map_point(r2_t *pP, r2x2_t *JP)
      {
        bool_t debug = debug_point(pP);
        
        /* Map native output image coords to user coords: */
        r2_map_projective(pP, &(osys.dir), pP, JP);
        if (debug) { r2_debug_point_jac("po", pP, JP, "\n"); }
        
        /* Apply the inverse radial correction: */
        if (o->kappa != 0) 
          { r2_t h = o->pixelSize;
            r2_map_radial(pP, &h, -(o->kappa), JP);
            if (debug) { r2_debug_point_jac("pr", pP, JP, "\n"); }
          }
        
        /* Map user coords to native input image coords: */
        r2_map_projective(pP, &(isys.inv), pP, JP);
        if (debug) { r2_debug_point_jac("pi", pP, JP, "\n"); }
        
        if (! o->extend)
          { /* Check domain: */
            bool_t invalid = FALSE;
            invalid |= ((pP->c[0] < 0) || (pP->c[0] > iCols));
            invalid |= ((pP->c[1] < 0) || (pP->c[1] > iRows));
            if (invalid) { pP->c[0] = pP->c[1] = NAN; }
          }
      }
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
    o->iCenter = TRUE;
    o->iOrg = (r2_t){{ 0.0, 0.0 }};
    o->oCenter = TRUE;
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
    
    argparser_get_keyword(pp, "-kappa");
    o->kappa = argparser_get_next_double(pp, -BIG, +BIG);
    
    argparser_get_keyword(pp, "-pixelSize");
    o->pixelSize.c[0] = argparser_get_next_double(pp, 1.0/BIG, BIG);
    o->pixelSize.c[1] = argparser_get_next_double(pp, 1.0/BIG, BIG);
    
    if (argparser_keyword_present(pp, "-interpolate"))
      { o->interpolate = (int32_t)argparser_get_next_int(pp, 0, 1); }
    else
      { /* The default depends on the input image: */
        o->interpolate = 0;
      }
    
    o->extend = argparser_keyword_present(pp, "-extend");
    
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (uint16_t)argparser_get_next_int(pp, 1, PNM_MAX_SAMPLE); }
    else
      { /* The default depends on the input image: */
        o->maxval = 0;
      }

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
    bool_t verbose
  )
  { uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    bool_t isMask = FALSE; /* Assume smooth distr of pixel values. */
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
    bool_t verbose
  )
  { int32_t chns = (int32_t)fim->sz[0];
    bool_t isMask = FALSE; /* Assume smooth distr of pixel values. */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, chns, NULL, NULL, NULL, maxval, TRUE, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    uint16_image_free(pim);
  }

