#define PROG_NAME "pnmprojmap"
#define PROG_DESC "applies a projective map to one or more pbm/ppm/pgm files"
#define PROG_VERS "2.0"

/* Last edited on 2023-10-15 17:07:49 by stolfi */

/* Copyright © 2002 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#include <pnmprojmap_info.h>
     
    /* !!! Add weights to features. !!! */
    /* !!! Add "-mapType" option !!! */
    /* !!! Add "-midway" option !!! */
    /* !!! Write out map matrix !!! */
    /* !!! Write out adjusted feature points !!! */

/* We need to set these in order to get {isnan}. What a crock... */
#undef __STRICT_ANSI__
#define _ISOC99_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <i2.h>
#include <r2.h>
#include <r2_extra.h>
#include <r2x2.h>
#include <r3x3.h>
#include <hr2.h>
#include <vec.h>
#include <fget.h>

#include <jspnm.h>
#include <jsstring.h>
#include <interval.h>
#include <ix.h>
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
#include <argparser_extra.h>
#include <argparser_geo.h>

#define USER_UNIT_MIN  (1.0e-100)
  /* A very small value, but still far from underflow. */
   
#define USER_UNIT_MAX  (1.0e+21)
  /* A very large value, but still far from overflow. */
 
#define USER_COORD_MAX  (1.0e+21)
  /* A very large value, but still far from overflow. */

#define IMG_SIZE_MAX (32*1024)
  /* A limit on image size, to avoid humongous mallocs. */

typedef struct feature_t
  { char *tag;   /* Key point tag. */
    r2_t pos;    /* Key point coordinates, in the image's input user coordinate system. */
    char *cmt;   /* The comment string for this point, or {NULL}. */
  } feature_t;
  /* Data about a key feature reaad from a ".pts" file. */

vec_typedef(feature_vec_t,feature_vec,feature_t);
  /* An expandable vector of {feature_t}. */

typedef struct image_options_t
  { char *name;           /* Input (and output) file name, minus prefix and extensions. */
    char *ext;            /* Input (and output) file extension. */
    char *matrix;         /* Name specified with the "map" attribute, or {NULL}. */
    char *points;         /* Name specified with the "points" attribute, or {NULL}. */
    /* User coordinates for input image: */
    double unit;          /* Input user coords unit, in pixels. */
    bool_t center;        /* If true, input origin is center; if false, use {org}. */
    r2_t org;             /* Input origin relative to default origin, if {!center}. */
  } image_options_t;

vec_typedef(image_options_vec_t,image_options_vec,image_options_t);

typedef struct options_t
  { /* User coordinate system axis directions: */
    bool_t yUp;                /* true' if the vertical axis points up, false otherwise. */
    bool_t xLeft;              /* true if the horizontal axis points left, false otherwise. */
    /* Input image data: */
    char *inPrefix;            /* Common prefix for all input file names. */
    image_options_vec_t image; /* The input image specs. */
    /* Map specification method: */
    bool_t fromPoints;         /* If true, the map is computed from matching points, else as matrix. */
    /* Output user coord system: */
    double oUnit;              /* Output user coords unit, in pixels. */
    bool_t oCenter;            /* If true, output origin is center; if false, use {oOrg}. */
    r2_t oOrg;                 /* Output origin relative to default origin, if {!oCenter}. */
    /* Output image size: */
    double oCols;              /* X size of output images (out user units), or -1 if not given. */
    double oRows;              /* Y size of output images (out user units), or -1 if not given. */
    /* Input image encoding and interpolation options: */
    bool_t isMask;              /* true to interpret samples as in masks. */
    bool_t extend;              /* True extends the the image by row/col replication. */
    float undef;                /* Input image epadding value, if {extend} is false. */
    int32_t interpolate;        /* Interpolation order. */
    /* Output image encoding: */
    bool_t noImages;            /* If true, skips computing and writing the output image. */
    char *outPrefix;            /* Common prefix for all output file names. */
    uint16_t maxval;            /* Output maxval (max int sample value) requested by user, or 0 if not given. */
    /* Debugging options: */
    bool_t verbose;             /* True to print global statistics. */
    i2_t debug;                 /* PBM indices of output pixel to debug; or (-1,-1). */
  } options_t;
  /* The transformation conceptually consists of moving the color from
    each point {ip} of the input image's domain (in the user's input
    coordinate system) to the point {op = M.dir(ip)} (in the user's
    output coordinate system). (In practice, we scan each output pixel
    {op} and set it to the input image's value in the neighborhood of
    point {ip = M.inv(op)}. */

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);
options_t *parse_options(int32_t argc, char **argv);
image_options_vec_t parse_image_options(argparser_t *pp);
void parse_image_center_org_unit_attributes(argparser_t *pp, image_options_t *im);

hr2_pmap_t compute_map_from_features(feature_vec_t *fta, feature_vec_t *ftb); 
  /* Finds all pairs of features {fa} in {fta} and {fb} in {ftb} that have
    the same tag. Then computes a projective map {M} that takes 
    {fa.pos} as close as possible to {fb.pos}, for every such pair. */

hr2_pmap_t best_map(feature_vec_t *fta, int32_vec_t *ixa, feature_vec_t *ftb, int32_vec_t *ixb);
  /* The index lists {ixa} and {ixb} must have the same length {nf} 
    and must contain the indices of features in {fta} and {ftb} with matching tags.
    The length {nf} should be more than the minimum number of points needed to determine
    the map. The procedure computes a map {N} such that the mean squared distance
    between {N(fta.e[ixa.e[k]])} and {N^{-1}(ftb.e[ixb.e[k]])} is as small as possible.
    Then returns the map {N^2}. */

void check_pmap(hr2_pmap_t *M, feature_vec_t *fta, int32_vec_t *ixa, feature_vec_t *ftb, int32_vec_t *ixb, double tol);
  /* The index lists {ixa} and {ixb} must be as in {best_map}. Checks
    the distances between {M(pa)} an {pb} for positions {pa,pb} of all pairs of features
    with matching tags. */

float_image_t *process_image
  ( options_t *o,
    image_options_t *oi,
    float_image_t *img,
    hr2_pmap_t *M_user
  );
  /* Applies the projective map {M} to the image {img}. */

void check_defines(void);
  /* Calls all the defined macros that are used in the {PROG_HELP} and {PROG_INFO} string,
    in bottom-up order, so that any bugs in them are easier to debug. */

void change_coord_systems
  ( hr2_pmap_t *P,    /* (IN/OUT) The projective transformation map. */
    hr2_pmap_t *isys, /* The map from pixel input to user input coordinates. */
    hr2_pmap_t *osys  /* The map from pixel output to user output coordinates. */
  );
  /* Assumes that {P} is a projective map from the input domain to the
    output domain, described in terms of the user-defined input and
    output coordinate systems. Modifies the map {P} so that it takes
    and returns points expressed in the pixel {float_image}
    coordinate systems, according to the given coordinate system options. */

feature_vec_t read_features(char *prefix, char *name);
  /* Reads the key points of one image from the file "{prefix}{name}.pts". */

hr2_pmap_t read_map_matrix(char *prefix, char *name);
  /* Reads a projective map as a 3x3 matrix (input user cooords to outout user coords)
     from the file "{prefix}{name}.mat". */

float_image_t *read_image
  ( char *prefix, 
    char *name,
    char *ext,
    int32_t *colsP,
    int32_t *rowsP,
    int32_t *chnsP,
    uint16_t *maxvalP,
    bool_t isMask,
    bool_t verbose
  );
  /* Reads a PBM/PGM/PPM image file from file "{prefix}{name}.{ext}", 
    converts it to a float image with samples in the range [0_1] and Y axis down. Returns the
    relevant image data. If {verbose} is true, prints image statistics
    to {stderr}. */
    
feature_vec_t read_features(char *prefix, char *name);
  /* Reads a list of feature points from file "{prefix}{name}.pts". */
        
int32_vec_t index_sort_features(feature_vec_t *ft);
  /* Let {nf} be {ft.ne}. Returns a vector of {nf} indices {ix} such
    that {ix.e[0..nf-1]} is a permutation of {0..nf-1]}, and the tags
    {ft.e[ix.e[0..nf-1]].tag} are in increasing lexical order. Fails
    with error if two distinct elements of {ft} have the same tag. */
    
void index_match_features(feature_vec_t *fta, int32_vec_t *ixa, feature_vec_t *ftb, int32_vec_t *ixb);
  /* Assumes that {ixa} and {ixb} are permutations of the indices of {fta} and {ftb}
    that indicate the order of their increasing tags, as returned by {index_sort_features}.
    Finds the {nf} elements of {fta} and {ftb} that have matching tags, and reduces
    both {ixa} and {ixb} to list only those elements, in the same order.  That is,
    after this call, {ixa.ne == ixb.ne == nf}, and {fta.e[ixa.e[k]].tag} and {ftb.e[ixb.e[k]].tag} are equal strings,
    for all {k} in {0..nf-1}. */
    
hr2_pmap_t read_map_matrix(char *prefix, char *name);
  /* Reads a projective map matrix from file "{prefix}{name}.map". */

void write_image
  ( char *prefix, 
    char *name,
    char *ext,
    float_image_t *fim,
    uint16_t maxval,
    bool_t isMask,
    bool_t verbose
  );
  /* Writes the float image {fim} to file "{prefix}{name}.{ext}" 
    as a PBM/PGM/PPM image file.
    Samples are converted fron the range [0_1]. Assumes Y axis down.
    Returns the relevant image data. If {verbose} is true, prints image
    statistics to {stderr}. */
    
char *make_filename(char *prefix, char *name, char *ext);
  /* Returns the string "{prefix}{name}.{ext}". */

/* IMPLEMENTATIONS */
    
void check_defines(void)
  { /* This procedure had better come before {main} so that #define bugs are caught upfront. */

    auto void ptc(char *name, char *t);
    
    char *tm1 = img_ops_HELP; ptc("tm1", tm1);
    char *t00 = imgc_parse_unit_HELP("-oUnit","_OUT"); ptc("t00", t00);
    char *t01 = imgc_parse_size_HELP("-oSize","_OUT"); ptc("t01", t01);
    char *t02 = PROG_HELP; ptc("t02", t02);
    char *t03 = PROG_INFO_COORDS_INTRO; ptc("t03", t03);
    char *t04 = PROG_INFO_DESC; ptc("t04", t04);
    char *t05 = xaxis_def; ptc("t05", t05);
    char *t06 = yaxis_def; ptc("t06", t06);
    char *t07 = yaxis_def; ptc("t07", t07);
    char *t08 = osize_def; ptc("t08", t08);
    char *t09 = PROG_INFO_MAP_FILES; ptc("t09", t09);
    char *t10 = PROG_INFO_POINTS_FILES; ptc("t10", t10);
    char *t11 = imgc_parse_x_axis_INFO_OPTS(xaxis_def); ptc("t11", t11);
    char *t12 = imgc_parse_y_axis_INFO_OPTS(yaxis_def); ptc("t12", t12);
    char *t13 = imgc_parse_unit_INFO_OPTS("unit","_IN[k]","each input image"); ptc("t13", t13);
    char *t14 = imgc_parse_center_org_INFO_OPTS("center","org","_IN[k]",org_def); ptc("t14", t14);
    char *t15 = imgc_unit_affects_org_INFO_OPTS("unit","org","the corresponding image"); ptc("t15", t15);
    char *t16 = imgc_parse_unit_INFO_OPTS("-oUnit","_OUT","every output image"); ptc("t16", t16);
    char *t17 = imgc_unit_affects_org_INFO_OPTS("-oUnit","org","the corresponding image"); ptc("t17", t17);
    char *t18 = imgc_unit_affects_org_INFO_OPTS("-oUnit","org","the corresponding image"); ptc("t18", t18);
    char *t19 = imgc_parse_size_INFO_OPTS("-oSize","_OUT","all output images",osize_def); ptc("t19", t19);
    char *t20 = PROG_INFO_OPTS; ptc("t20", t20);
    char *t21 = PROG_INFO; ptc("t21", t21);
    
    void ptc(char *name, char *text)
      { fprintf(stderr, "--- %s ----------------------------------------------------------------\n", name);
        argparser_print_text(stderr, text, 72);
        fprintf(stderr, "------------------------------------------------------------------------\n");
      }
  }
    
float_image_t *process_image
  ( options_t *o,
    image_options_t *oi,
    float_image_t *img_in,
    hr2_pmap_t *M_user
  )
  {
    /* Get input image channels and size (in pixels): */
    int32_t chns, iCols, iRows;
    float_image_get_size(img_in, &chns, &iCols, &iRows);

    /* Compute output image size in pixels: */
    int32_t oCols, oRows;
    imgc_compute_output_size_in_pixels(iCols, iRows, oi->unit, o->oCols, o->oRows, o->oUnit, &oCols, &oRows, IMG_SIZE_MAX);
    if (o->verbose) { fprintf(stderr, "output image size (pixels) = %d %d\n", oCols, oRows); }

    /* Compute pixel-to-user coord system map {isys} for input image: */
    hr2_pmap_t isys = imgc_coord_sys_map(o->xLeft, o->yUp, oi->unit, oi->center, &(oi->org), iCols, iRows);
    if (o->verbose) { imgc_print_pmap(stderr, "input pixel", "input user", "isys", &(isys)); }

    /* Compute pixel-to-user coord system map {osys} for the output image: */
    hr2_pmap_t osys = imgc_coord_sys_map(o->xLeft, o->yUp, o->oUnit, o->oCenter, &(o->oOrg), oCols, oRows);
    if (o->verbose) { imgc_print_pmap(stderr, "output pixel", "output user", "osys", &(osys)); }

    /* Build the pixel-to-pixel projective map {M_pix}: */
    hr2_pmap_t M_pix = (*M_user);
    change_coord_systems(&M_pix, &isys, &osys);
    if (o->verbose) { imgc_print_pmap(stderr, "input pixel", "output pixel", "M_pix", &M_pix); }
    /* Check whether tha projective map is the identity: */
    bool_t ident_proj = hr2_pmap_is_identity(&M_pix);
 
    bool_t debugging = ((o->debug.c[0] >= 0) && (o->debug.c[1] >= 0));
    if (debugging && o->verbose)
      { /* Map center of debugging pixel from pixel coords to pixel output coords {odb_pix}: */
        r2_t odb_pix, odb_usr;
        odb_pix.c[0] = o->debug.c[0] + 0.5;
        odb_pix.c[1] = o->debug.c[1] + 0.5;
        r2_map_projective(&odb_pix, &(osys.dir), &odb_usr, NULL);
        /* Map center of debugging pixel from pixel coords to pixel input coords {idb_pix}: */
        r2_t idb_usr, idb_pix;
        r2_map_projective(&odb_usr, &(M_user->inv), &idb_usr, NULL);
        r2_map_projective(&idb_usr, &(isys.inv), &idb_pix, NULL);
        i2_t idb_pix_ix = (i2_t){{ (int32_t)floor(idb_pix.c[0]), (int32_t)floor(idb_pix.c[1]) }};
        fprintf(stderr, "watching output pixel col = %d row = %d whose center is\n", o->debug.c[0], o->debug.c[1]);
        fprintf(stderr, "  = (%7.1f,%7.1f) pixel in output pixel coords\n", odb_pix.c[0], odb_pix.c[1]);
        fprintf(stderr, "  = (%7.1f,%7.1f) pixel in output user coords\n", odb_usr.c[0], odb_usr.c[1]);
        fprintf(stderr, "  = (%7.1f,%7.1f) pixel in input user coords\n", idb_usr.c[0], idb_usr.c[1]);
        fprintf(stderr, "  = (%7.1f,%7.1f) pixel in input pixel coords\n", idb_pix.c[0], idb_pix.c[1]);
        fprintf(stderr, "  inside input pixel col = %d row = %d\n", idb_pix_ix.c[0], idb_pix_ix.c[1]);
      }

    /* Allocate output image: */
    float_image_t *img_ot = float_image_new(chns, oCols, oRows);

    /* Apply the transformation from input image to output image: */

    auto void map_point(r2_t *pP, r2x2_t *JP);
      /* Maps a point {*p} of the OUTPUT image (in the {float_image}
        pixel coordinates, with (0,0) at bottom left) to the
        corresponding point of the INPUT image and stores that into {*p}. If the point
        falls outside the input image's domain, and {o->extend} is
        false, sets {*p} to {(NAN,NAN)}.  Also multiplies {*J} by the
        Jacobian od the map. */

    auto bool_t debug_point(r2_t *pP);
      /* True if {*pP} is the watched pixel. */
      
    if (o->verbose) { fprintf(stderr, "transforming the image ...\n"); }
    ix_reduction_t red = ix_reduction_SINGLE;
    float_image_transform_all(img_in, red, &map_point, o->undef, TRUE, o->interpolate, &debug_point, img_ot);

    return img_ot;

    bool_t debug_point(r2_t *pP)
      { return debugging &&
          (fabs(pP->c[0] + 0.5 - o->debug.c[0]) <= 0.499) &&
          (fabs(pP->c[1] + 0.5 - o->debug.c[1]) <= 0.499);
      }

    void map_point(r2_t *pP, r2x2_t *JP)
      {
        bool_t debug = debug_point(pP);

        /* Apply the inverse projective map to the ouput image point: */
        if (! ident_proj)
          { r2_map_projective(pP, &(M_pix.inv), pP, JP);
            if (debug) { r2_debug_point_jac("pp", pP, JP, "\n"); }
          }

        if (! o->extend)
          { /* Check domain: */
            bool_t invalid = FALSE;
            invalid |= ((pP->c[0] < 0) || (pP->c[0] > iCols));
            invalid |= ((pP->c[1] < 0) || (pP->c[1] > iRows));
            if (invalid) { pP->c[0] = pP->c[1] = NAN; }
          }
      }
  }

int32_t main(int32_t argc, char **argv)
  {
    /* check_defines(); */

    /* Parse command line options: */
    options_t *o = parse_options(argc, argv);
    
    feature_vec_t ft_ref = feature_vec_new(0); /* Reference point set, if any. */
    int32_t ni = o->image.ne; /* Number of input images. */
    for (int32_t ki = 0; ki < ni; ki++)
      { image_options_t *oik = &(o->image.e[ki]);
        hr2_pmap_t Mk_user;
        if (o->fromPoints)
          { char *pts_name = (oik->points != NULL ? oik->points : oik->name);
            feature_vec_t ftk = read_features(o->inPrefix, pts_name); 
            if (ki == 0)
              { ft_ref = ftk; Mk_user = hr2_pmap_identity(); }
            else
              { Mk_user = compute_map_from_features(&ftk, &ft_ref); }
          }
        else
          { char *map_name = (oik->matrix != NULL ? oik->matrix : oik->name);
            Mk_user = read_map_matrix(o->inPrefix, map_name);
          }
        if (o->verbose) { imgc_print_pmap(stderr, "input user", "output user", "M", &(Mk_user)); }

        if (! o->noImages)
          { int32_t iCols, iRows, chns;
            uint16_t maxval_in;
            float_image_t *fim_in = read_image
              ( o->inPrefix, oik->name, oik->ext, 
                &iCols, &iRows, &chns, &maxval_in, 
                o->isMask, o->verbose
              );
            uint16_t maxval_out = o->maxval;
            if (maxval_out == 0)
              { uint16_t mvmin = (chns == 3 ? 255 : PNM_MAX_SAMPLE); 
                maxval_out = (maxval_in < mvmin ? mvmin : maxval_in);
              }
            float_image_t *fim_out = process_image(o, oik, fim_in, &Mk_user);
            write_image
              ( o->outPrefix, oik->name, oik->ext,
                fim_out, maxval_out, o->isMask, o->verbose
              ); 
          }
      }
 
    if (o->verbose) { fprintf(stderr, "done.\n"); }
    return 0;
  }
 
char *make_filename(char *prefix, char *name, char *ext)
  { char *fname = NULL;
    asprintf(&fname, "%s%s.%s", prefix, name, ext);
    return fname;
  }

void change_coord_systems
  ( hr2_pmap_t *P,
    hr2_pmap_t *isys, /* The map from pixel input to user input coordinates. */
    hr2_pmap_t *osys  /* The map from pixel output to user output coordinates. */
  )
  {
    /* Replace {P} by {(isys)*P*(osys^-1): */
    (*P) = hr2_pmap_compose(isys, P);
    hr2_pmap_t syso = hr2_pmap_inv(osys);
    (*P) = hr2_pmap_compose(P, &syso);
  }

float_image_t *read_image
  ( char *prefix,
    char *name, 
    char *ext,
    int32_t *colsP,
    int32_t *rowsP,
    int32_t *chnsP,
    uint16_t *maxvalP,
    bool_t isMask,
    bool_t verbose
  )
  { char *fname = make_filename(prefix, name, ext);
    FILE *rd = open_read(fname, verbose);
    free(fname);
    uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    fclose(rd);
    bool_t yup = FALSE; /* Use the traditional image system to reduce confusion. */
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, yup, verbose);
    (*colsP) = pim->cols;
    (*rowsP) = pim->rows;
    (*chnsP) = pim->chns;
    (*maxvalP) = pim->maxval;
    uint16_image_free(pim);
    return fim;
  }

feature_vec_t read_features(char *prefix, char *name)
  { char *fname = make_filename(prefix, name, "pts");
    FILE *rd = open_read(fname, TRUE);
    feature_vec_t ft = feature_vec_new(0);
    int32_t nf = 0; /* Number of points read. */
    int32_t line_num = 0; /* Number of lines found. */

    auto void mp_error(char *tag, char *msg);
    auto r2_t read_r2(char *tag);
    
    while (TRUE)
      { fget_skip_spaces(rd);
        if (fget_test_eof(rd)) { break; }
        line_num++; 
        if (fget_test_comment_or_eol(rd, '#', NULL)) { continue; }
        char *tag = fget_string(rd);
        uint8_t ch = (uint8_t)tag[0];
        if (((ch < 'a') || (ch > 'z')) && ((ch < 'A') || (ch > 'Z')))
          { mp_error(tag, " tag should start with letter"); }
        r2_t pos = read_r2(tag);
        fget_skip_spaces(rd);
        if (fget_test_char(rd, 'C'))
          { /* Skip feature radius and angle: */
            /* double rad = */ (void)fget_double(rd);
            /* double ang = */ (void)fget_double(rd);
          }
        else if (fget_test_char(rd, 'E'))
          { /* Skip conjugated vectors of feature ellipsoid: */
            /* r2_t u = */ (void)read_r2(tag);
            /* r2_t v = */ (void)read_r2(tag);
          }
        else if (fget_test_char(rd, 'P'))
          { /* No more data */
          }
        else 
          { /* Should not happen: */
            mp_error(tag, " no feature type letter");
          }
        fget_skip_spaces(rd);
        char *cmt = NULL;
        fget_comment_or_eol(rd, '#', &cmt);
        if (cmt != NULL)
          { char *t = trim_spaces(cmt, TRUE, TRUE);
            free(cmt); 
            cmt = t;
          }
        if ((cmt != NULL) && (strlen(cmt) == 0)) { cmt = NULL; }
            
        feature_vec_expand(&ft, nf);
        feature_t *ftk = &(ft.e[nf]);
        (*ftk) = (feature_t){ .tag = tag, .pos = pos, .cmt = cmt };
        nf++;
      }
    feature_vec_trim(&ft, nf);
    free(fname);
    fclose(rd);
    return ft;
    
    void mp_error(char *tag, char *msg)
      { fprintf(stderr, "%s:%d: **", fname, line_num);
        if (tag != NULL) { fprintf(stderr, " \"%s\":", tag); }
        fprintf(stderr, " %s\n", msg);
        affirm(FALSE, "aborted");
      }

    r2_t read_r2(char *tag)
      { r2_t p;
        for (int32_t j = 0; j < 2; j++)
          { p.c[j] = fget_double(rd);
             if (fabs(p.c[j]) > USER_COORD_MAX) { mp_error(tag, " coord too big"); }
          }
        return p;
      }
  }
       
hr2_pmap_t read_map_matrix(char *prefix, char *name)
  { char *fname = make_filename(prefix, name, "map");
    FILE *rd = open_read(fname, TRUE);
    free(fname);
    hr2_pmap_t M;
    int32_t line_num = 0;
    /* Skip blank and comment lines: */
    while (TRUE)
      { fget_skip_spaces(rd);
        if (fget_test_eof(rd)) { break; }
        line_num++; 
        fprintf(stderr, "  (a) %d\n", line_num);
        if (! fget_test_comment_or_eol(rd, '#', NULL)) { break; }
      }
    /* Read the direct matrix in three lines: */
    for (int32_t i = 0; i < 3; i++)
      { fprintf(stderr, "  (b) %d\n", line_num);
        for (int32_t j = 0; j < 3; j++)
          { M.dir.c[i][j] = fget_double(rd); }
        demand(fget_test_comment_or_eol(rd, '#', NULL), "spurious data at end of line");
        line_num++;
      }
    r3x3_inv(&(M.dir), &(M.inv));
    return M;
  }

void write_image
  ( char *prefix,
    char *name, 
    char *ext,
    float_image_t *fim,
    uint16_t maxval,
    bool_t isMask,
    bool_t verbose
  )
  { int32_t chns = (int32_t)fim->sz[0];
    bool_t yup = FALSE; /* Use the traditional image system to for consistency with {read_image}. */
    uint16_image_t *pim = float_image_to_uint16_image(fim, isMask, chns, NULL, NULL, NULL, maxval, yup, verbose);
    bool_t forceplain = FALSE;
    char *fname = make_filename(prefix, name, ext);
    FILE *wr = open_write(fname, verbose);
    free(fname);
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose); 
    fclose(wr);
    uint16_image_free(pim);
  }

hr2_pmap_t compute_map_from_features(feature_vec_t *fta, feature_vec_t *ftb)
  { 
    int32_vec_t ixa = index_sort_features(fta);
    int32_vec_t ixb = index_sort_features(ftb);
    index_match_features(fta, &ixa, ftb, &ixb);
    assert(ixa.ne == ixb.ne);
    int32_t nf = ixa.ne; /* Number of matched features. */
    
    /* Extract the coordinate of the matched features: */
    r2_t *pa = talloc(nf, r2_t);
    r2_t *pb = talloc(nf, r2_t);
    for (int32_t kf = 0; kf < nf; kf++)
      { pa[kf] = fta->e[ixa.e[kf]].pos;
        pb[kf] = ftb->e[ixb.e[kf]].pos;
      }
        from many points
        
    int32_t nf_min = 4; /* Min points for "exact" map computation. */
    hr2_pmap_t M;
    if (nf > nf_min)
      { /* Too many features for simple formula, must use optimization: */
        M = best_map(fta, &ixa, ftb, &ixb);
      }
    else if (nf == 0)
      { M = hr2_pmap_identity(); }
    else if (nf == 1)
      { /* Translation */ 
        r2_t pa = fta->e[ixa.e[0]].pos;
        r2_t pb = ftb->e[ixb.e[0]].pos;
        r2_t d; r2_sub(&pb, &pa, &d);
        M = hr2_pmap_translation(&d);
      }
    else
      { hr2_pmap_t MA, MB; /* The two halves of the map. */
        if (nf == 2)
          { /* Similarity */
            r2_t pa = fta->e[ixa.e[0]].pos;
            r2_t qa = fta->e[ixa.e[1]].pos;
            MA = hr2_pmap_similarity_from_two_points(&pa, &qa, FALSE);

            r2_t pb = ftb->e[ixb.e[0]].pos;
            r2_t qb = ftb->e[ixb.e[1]].pos;
            MB = hr2_pmap_similarity_from_two_points(&pb, &qb, FALSE);
          }
        else if (nf == 3)
          { /* Affine */
            r2_t pa = fta->e[ixa.e[0]].pos;
            r2_t qa = fta->e[ixa.e[1]].pos;
            r2_t ra = fta->e[ixa.e[2]].pos;
            MA = hr2_pmap_aff_from_three_points(&pa, &qa, &ra);

            r2_t pb = ftb->e[ixb.e[0]].pos;
            r2_t qb = ftb->e[ixb.e[1]].pos;
            r2_t rb = ftb->e[ixb.e[2]].pos;
            MB = hr2_pmap_aff_from_three_points(&pb, &qb, &rb);
          }
        else if (nf == 4)
          { /* General projective. */
            hr2_point_t pa = hr2_from_r2(&(fta->e[ixa.e[0]].pos));
            hr2_point_t qa = hr2_from_r2(&(fta->e[ixa.e[1]].pos));
            hr2_point_t ra = hr2_from_r2(&(fta->e[ixa.e[2]].pos));
            hr2_point_t sa = hr2_from_r2(&(fta->e[ixa.e[3]].pos));
            MA = hr2_pmap_from_four_points(&pa, &qa, &ra, &sa);

            hr2_point_t pb = hr2_from_r2(&(ftb->e[ixb.e[0]].pos));
            hr2_point_t qb = hr2_from_r2(&(ftb->e[ixb.e[1]].pos));
            hr2_point_t rb = hr2_from_r2(&(ftb->e[ixb.e[2]].pos));
            hr2_point_t sb = hr2_from_r2(&(ftb->e[ixb.e[3]].pos));
            MB = hr2_pmap_from_four_points(&pb, &qb, &rb, &sb);
          }
        else
          { assert(FALSE); }
        M = hr2_pmap_inv_compose(&MA, &MB);
      }
    
    /* !!! The tolerance {tol} should be a user parameter. */
    check_pmap(&M, fta, &ixa, ftb, &ixb, 1.0);
    return M;
  }

hr2_pmap_t best_map(feature_vec_t *fta, int32_vec_t *ixa, feature_vec_t *ftb, int32_vec_t *ixb)
  { demand(ixa->ne == ixb->ne, "counts of matched points differ");
    int32_t nf = ixa->ne; /* Number of matched features. */
    demand(nf <= fta->ne, "bad {fta.ne}");
    demand(nf <= ftb->ne, "bad {ftb.ne}");
    
    hr2_pmap_t M = initial_map(fta, ixa, ftb, ixb);
    
    demand(FALSE, "NOT IMPLEMENTED");
    return hr2_pmap_identity();
  }
    
hr2_pmap_t initial_map(feature_vec_t *fta, int32_vec_t *ixa, feature_vec_t *ftb, int32_vec_t *ixb)
  { demand(ixa->ne == ixb->ne, "counts of matched points differ");
    int32_t nf = ixa->ne; /* Number of matched features. */
    demand(nf <= fta->ne, "bad {fta.ne}");
    demand(nf <= ftb->ne, "bad {ftb.ne}");
    
    /* Pick a minimum subset of key points, trying to get them as spread out as possible: */
    int32_t nf_min = 4; /* Min number of points for driect solution. */
    demand(nf >= nf_min, "{nf} too small");
    r2_t *pa[nf_min], *pb[nf_min]; /* Selected points for initial map. */
    select_some_features(fta, ixa, ftb, ixb, nf_min, pa, pb);
    /* Compute the map from those points: */
    assert(nf_min == 4);
    hr2_pmap_t Ma = hr2_pmap_from_four_points(pa[0], pa[1], pa[2], pa[3]);
    hr2_pmap_t Mb = hr2_pmap_from_four_points(pb[0], pb[1], pb[2], pb[3]);
    hr2_pmap_t M = hr2_pmap_inv_compose(&Ma, &Mb);
    return M;
  }
    
void select_some_features
  ( int32_t nf,
    r2_t pa[],
    r2_t pb[],
    int32_t ns,
    int32_t kf_sel[]
  )
  { demand(ns >= 0, "invalid {ns}");
    demand(ns <= nf, "not enough matched features");
    
    if (ns == 0) 
      { return; }
    else if (ns == nf) 
      { /* Must take all of them: */
        for (int32_t ks = 0; ks < ns; ks++)
          { pa[ks] = &(fta.e[ixa.e[ks]].pos);
            pb[ks] = &(fta.e[ixb.e[ks]].pos);
          }
      }
    else 
      { r2_t ca, cb; /* Barycenters of the two feature sets. */
        double sza, szb; /* Estimated scale lengths of the two sets. */
        compute_feature_center_scale(fta, ixa, &ca, &sza);
        compute_feature_center_scale(ftb, ixb, &cb, &szb);
        if (ns == 1)
          { /* Pick the feature that on both sets seems to be closer to the center: */
            int32_t k0 = select_extremal_feature(fta, ixa, &ca, sza, ftb, ixb, &cb, szb, -1);
            pa[0] = &(fta.e[ixa.e[k0]].pos);
            pb[0] = &(fta.e[ixb.e[k0]].pos);
          }
        else
          { /* Start with the most outlying feature: */ 
            int32_t k0 = select_extremal_feature(fta, ixa, &ca, sza, ftb, ixb, &cb, szb, +1);
            pa[0] = &(fta.e[ixa.e[k0]].pos);
            pb[0] = &(fta.e[ixb.e[k0]].pos);
            /* Pick additional features that are as far as possible from each other: */
            for (int32_t js = 1; js < ns; js++)
              { int32_t kj = select_most_isolated_feature(js, fta, ixa, pa, sza, ftb, ixb, pb, szb);
                pa[js] = &(fta.e[ixa.e[kj]].pos);
                pb[js] = &(ftb.e[ixb.e[kj]].pos);
              }
          }
      }
  }
        
void compute_feature_center_scale(int32_t nf, r2_t p[], r2_t *c_P, double *sz_P)
  { 
    demand(nf >= 2, "needs at least two points");
    
    r2_t sum_p = (r2_t){{ 0, 0 }}; /* Sum of all matched feature positions. */
    double sum_d2 = 0; /* Sum of squares of distances of pairs of matched freatures. */
    double
    for (int32_t kf = 0; kf < nf; kf++)
      { r2_t *pk = &(p[kf]);
        r2_add(pk, &sum_p, &sum_p);
        for (int32_t jf = 0; jf < kf; jf++)
          { r2_t *pj = &(p[jf]);
            double d2 = r2_dist_sqr(pk, pj);
            sum_d2 += d2;
          }
      }
    /* Compute barycenter: */
    r2_scale(1.0/nf, &sum_pos, c_P);
    /* Compute estimated scale: */
    (*sz_P) = sqrt(sum_d2/(nf*(nf-1)));

  }
  
int32_t select_extremal_feature
  ( int32_t nf,
    r2_t pa[],
    r2_t *ca,
    double sza,
    r2_t pb[],
    r2_t *cb,
    double szb,
    sign_t dir
  )
  { demand(nf >= 1, "needs at least one feature");
    
    if (nf == 2)
      { /* Either one is the centermost one. */
        return 0;
      }
    else
      { double rdmax = -INF; /* Min relative distance from center, times {sgn}. */
        int32_t kfmax = -1;
        for (int32_t kf = 0; kf < nf; kf++)
          { double rda = r2_dist(&(pa[kf]), ca)/sza;
            double rdb = r2_dist(&(pb[kf]), cb)/szb;
            doubel rd = sgn*(rda + rdb)/2;
            if (rd > rdmax) { kfmax = kf; rdmax = rd; }
          }
        assert(kfmax >= 0);
        return kfmax;
      }
  }
  
int32_t select_most_isolated_feature
  ( int32_t ns,
    int32_t kf_sel[],
    int32_t nf
    r2_t pta[],
    double sza,
    r2_t ptb[],
    double szb
  )
  { demand(ns >= 1, "furthest from what?");
    demand(ns < nf, "No more features");
    
    auto void double min_rel_dist(int32_t kf, r2_t pt[], double double sz);
      /* Returns the min distance between the matched feature {pt[kf]]} and each of
        the already selected matched features {pt[jf]} where {jf} is any of {kf_sel[0..ns-1]},
        scaled by {1/sz}.  In particular, returns 0 if {kf} is one of those already selected features. */ 
    
    double rdmax = -INF; /* Max relative mean distance from previus ones. */
    int32_t kfmax = -1;
    for (int32_t kf = 0; kf < nf; kf++)
      { double rda = min_rel_dist(kf, pta, sza);
        double rdb = min_rel_dist(kf, ptb, szb);
        double rd = rda + rdb)/2;
        if (rd > rdmax) { kfmax = kf; rdmax = rd; }
      }
    demand(kfmax >= 0, "can't find another distinct point");
    return kfmax;

    double min_rel_dist(int32_t kf, r2_t pt[], double double sz)
      { double rdmin = +INF;
        for (int32_t js = 0; js < ns; js++)
          { int32_t jf = kf_sel[js];
            if (jf == kf) { /* Feature {kf} aready selected: */ return 0.0; }
            double rd = r2_dist(&(pt[kf]), &(pt[jf]))/sz;
            if (rd == 0) { return 0; }
            if (rd < rdmin) { rdmin = rd; }
          }
        return rdmin;
      }
  }

void check_pmap(hr2_pmap_t *M, feature_vec_t *fta, int32_vec_t *ixa, feature_vec_t *ftb, int32_vec_t *ixb, double tol)
  { demand(ixa->ne == ixb->ne, "mismatched lengths {ixa,ixb}"); /* Paranoia. */
    int32_t nf = ixa->ne;
    
    auto void prbug(char *name0, r2_t *p0, char *dir, r2_t *q0, char *name1, r2_t *p1, double d);
      /* Prints to {stderr} the data for a pair whose mapping error exceeds the tolerance. */
      
    auto void prtp(r2_t *p);
      /* Prints  to {stderr} the point {p} in a fixed-width format. */
    
    for (int32_t k = 0; k < nf; k++)
      { int32_t ka = ixa->e[k];
        demand((ka >= 0) && (ka < fta->ne), "bad {ixa} index vector"); /* Paranoia. */
        char *taga = fta->e[ka].tag;
        r2_t *pa = &(fta->e[ka].pos);
        
        int32_t kb = ixb->e[k];
        demand((kb >= 0) && (kb < ftb->ne), "bad {ixb} index vector"); /* Paranoia. */
        char *tagb = ftb->e[kb].tag;
        r2_t *pb = &(ftb->e[kb].pos);
      
        demand(strcmp(taga, tagb) == 0, "mismatched tags");  /* Paranoia. */
        
        /* !!! Should use midway distance !!! */
        
        r2_t qa = hr2_pmap_r2_point(pa, M);
        double d_dir = r2_dist(&qa, pb);

        r2_t qb = hr2_pmap_inv_r2_point(pb, M);
        double d_inv = r2_dist(&qb, pa);
        
        if ((d_dir > tol) || (d_inv > tol))
          { fprintf(stderr, "!!  mapping error too large for points \"%s\"\n", taga);
            prbug("a", pa, "dir", &qa, "b", pb, d_dir);
            prbug("p", pb, "inv", &qb, "a", pa, d_inv);
          }
      }

    return;
          
    void prbug(char *name0, r2_t *p0, char *dir, r2_t *q0, char *name1, r2_t *p1, double d)
      { fprintf(stderr, "    p%s = ", name0); prtp(p0); 
        fprintf(stderr, " -%s-> ", dir); prtp(q0); 
        fprintf(stderr, " p%s = ", name1); prtp(p1);
        fprintf(stderr, " dist = %.2f\n", d);
      }
      
    void prtp(r2_t *p)
      { r2_gen_print(stderr, p, "%7.2f", "( ", " ", " )"); }
        
  }
        

int32_vec_t index_sort_features(feature_vec_t *ft)
  { 
    int32_vec_t ix = int32_vec_new(ft->ne);
    for (int32_t i = 0; i < ft->ne; i++) { ix.e[i] = i; }
    
    auto int32_t cmptag(const void *p, const void *q);
    
    qsort((void *)ix.e, (size_t)ix.ne, sizeof(int32_t), &cmptag);
    return ix;
    
    int32_t cmptag(const void *p, const void *q)
      { int32_t ip = *((int32_t*)p); char *tagp = ft->e[ip].tag;
        int32_t iq = *((int32_t*)q); char *tagq = ft->e[iq].tag;
        int32_t cmp = strcmp(tagp, tagq);
        demand(cmp != 0, "duplicated feature tag");
        return cmp;
      }
  }

void index_match_features(feature_vec_t *fta, int32_vec_t *ixa, feature_vec_t *ftb, int32_vec_t *ixb)
  { int32_t nfa = fta->ne; demand(ixa->ne == nfa, "wrong index vector size");
    int32_t nfb = ftb->ne; demand(ixb->ne == nfb, "wrong index vector size");
    int32_t nf = 0; /* Number of matched features. */
    int32_t ka = 0;
    int32_t kb = 0;
    while ((ka < nfa) || (kb < nfb))
      { int32_t cmp;
        if (ka > nfa)
          { cmp = +1; }
        else if (kb > nfb) 
          { cmp = -1; }
        else
          { assert((ka < nfa) && (kb < nfb));
            char *taga = fta->e[ixa->e[ka]].tag;
            char *tagb = fta->e[ixa->e[ka]].tag;
            cmp = strcmp(taga, tagb);
          }
        if (cmp < 0)
          { /* Discard unmatched feature from {fta}: */ assert(ka < nfa); ka++; }
        else if (cmp > 0)
          { /* Discard unmatched feature from {ftb}: */ assert(kb < nfb); kb++; }
        else
          { /* Keep matched feature: */ 
            assert((nf <= ka) && (ka < nfa));
            assert((nf <= kb) && (kb < nfb));
            ixa->e[nf] = ixa->e[ka]; ka++; 
            ixb->e[nf] = ixb->e[kb]; kb++; 
            nf++;
          }
      }
    int32_vec_trim(ixa, nf);
    int32_vec_trim(ixb, nf);
  }

options_t *parse_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);

    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");

    /* Set defaults for input and output coord systems: */
    o->xLeft = FALSE;
    imgc_parse_x_axis(pp, &(o->xLeft));
    
    o->yUp = FALSE;
    imgc_parse_y_axis(pp, &(o->yUp));
   
    argparser_get_keyword(pp, "-inPrefix");
    o->inPrefix = argparser_get_next_non_keyword(pp);

    o->image = parse_image_options(pp);
    
    /* Parse output user coord system specs: */
    imgc_parse_unit(pp, "-oUnit", &(o->oUnit));
    imgc_parse_center_org(pp, "-oCenter", &(o->oCenter), "-oOrg", &(o->oOrg));

    /* The default output size depends on the input image size, so leave {-1}: */
    o->oCols = -1.0;
    o->oRows = -1.0;
    imgc_parse_size(pp, "-oSize", &(o->oCols), &(o->oRows));

    if (argparser_keyword_present(pp, "-fromMatrix"))
      { o->fromPoints = FALSE; }
    else if (argparser_keyword_present(pp, "-fromPoints"))
      { o->fromPoints = TRUE; }
    else
      { argparser_error(pp, "must specify \"-fromMatrix\" or \"-fromPoints\""); }
      
    o->noImages = argparser_keyword_present(pp, "-noImages");
   
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    if (argparser_keyword_present(pp, "-interpolate"))
      { o->interpolate = (int32_t)argparser_get_next_int(pp, 0, 1); }
    else
      { o->interpolate = 0; }

    if (argparser_keyword_present(pp, "-extend"))
      { o->extend = TRUE; o->undef = 0.5f; }
    else if (argparser_keyword_present(pp, "-undef"))
      { o->extend = FALSE;
        o->undef = (float)argparser_get_next_double(pp, 0.0, 1.0);
      }
    else
      { o->extend = FALSE; o->undef = 0.5f; }

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
      { o->debug.c[0] = (int32_t)argparser_get_next_int(pp, -1, IMG_SIZE_MAX);
        o->debug.c[1] = (int32_t)argparser_get_next_int(pp, -1, IMG_SIZE_MAX);
      }
    else
      { o->debug = (i2_t){{ -1, -1 }}; }

    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }

image_options_vec_t parse_image_options(argparser_t *pp)
  { 
    image_options_vec_t image = image_options_vec_new(10);
    
    int32_t ni = 0; /* Number of "-image" options given. */
    while (argparser_keyword_present(pp, "-image"))
      { image_options_vec_expand((&image), ni);
        image_options_t *oi = &(image.e[ni]);
        oi->name = argparser_get_next_non_keyword(pp);
        oi->ext = argparser_get_next_non_keyword(pp);
        parse_image_center_org_unit_attributes(pp, oi);
        ni++;
      }
    image_options_vec_trim((&image), ni);
    return image;
  }

void parse_image_center_org_unit_attributes(argparser_t *pp, image_options_t *oi)
  { oi->center = FALSE;
    oi->org = (r2_t){{ 0.0, 0.0 }};
    oi->unit = 1.0;
    oi->matrix = NULL;
    oi->points = NULL;
    bool_t center_org_given = FALSE; /* "center" or "org" specified for this image. */
    bool_t unit_given = FALSE; /* "unit" specified for this image. */
    while (argparser_next_is_non_keyword(pp))
      { char *subop = argparser_get_next_non_keyword(pp);
        if (strcmp(subop, "center") == 0)
          { if (center_org_given) 
              { argparser_error(pp, "repeated or contradictory \"center\"/\"org\" attributes"); }
            oi->center = TRUE;
            center_org_given = TRUE;
          }
        else if (strcmp(subop, "org") == 0)
          { if (center_org_given) 
              { argparser_error(pp, "repeated or contradictory \"center\"/\"org\" attributes"); }
            oi->org.c[0] = argparser_get_next_double(pp, -USER_COORD_MAX, +USER_COORD_MAX);
            oi->org.c[1] = argparser_get_next_double(pp, -USER_COORD_MAX, +USER_COORD_MAX);
            center_org_given = TRUE;
          }
        else if (strcmp(subop, "unit") == 0)
          { if (unit_given) { argparser_error(pp, "repeated \"unit\" attribute"); }
            oi->unit = argparser_get_next_double(pp, USER_UNIT_MIN, USER_UNIT_MAX);
            unit_given = TRUE;
          }
        else if (strcmp(subop, "matrix") == 0)
          { if (oi->matrix != NULL) { argparser_error(pp, "repeated \"matrix\" attribute"); }
            oi->matrix = argparser_get_next_non_keyword(pp);
          }
        else if (strcmp(subop, "points") == 0)
          { if (oi->points != NULL) { argparser_error(pp, "repeated \"points\" attribute"); }
            oi->points = argparser_get_next_non_keyword(pp);
          }
        else
          { argparser_error(pp, "invalid image attribute"); }
      }
  }

vec_typeimpl(feature_vec_t,feature_vec,feature_t);
vec_typeimpl(image_options_vec_t,image_options_vec,image_options_t);
