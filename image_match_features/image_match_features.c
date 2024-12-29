#define PROG_NAME "image_match_features"
#define PROG_DESC "precisely locates matching image features in two images"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-21 14:00:54 by stolfi */

#define image_match_features_C_COPYRIGHT \
    "© 2020 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -image1 {FILE1} {LX1} {LY1} {HX1} {HY1} \\\n" \
  "    -image2 {FILE1} {LX2} {LY2} {HX2} {HY2} \\\n" \
  "    [ -fix {IFIX} ] \\\n" \
  "    [ -maxIter {MAX_ITER} ] \\\n" \
  "    [ -maxAdjust {MAX_ADJ} ] \\\n" \
  "    [ -writeFeatures {OUT_PREFIX} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    < {IN_PAIR_FILE} \\\n" \
  "    > {OUT_PAIR_FILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "   Reads from standard input a file of approximately corresponding" \
  " points in two given images.  Independently refines the positions of each" \
  " point pair so as to maximize the similarity of their" \
  " neighborhoods.  Writes the refined positions to standard output.\n" \
  "\n" \
  "   The points had better be coordinates of well-matcheable features, that" \
  " is, details that become quite different when translated by any" \
  " amount.  Thus spots, corners, crosses, or butterfly patterns are" \
  " better than long straight edges or lines.\n" \
  "\n" \
  "  For convenience, all feature positions are specified in arbitrary" \
  " Client {X,Y} coordinates, instead of image pixel column and row" \
  " indices.  This makes it easier to run the program on variously" \
  " scaled versions of the images.\n" \
  "\n" \
  "INPUT FILE.\n" \
  "  " imft_data_format_INFO "\n" \
  "\n" \
  "OUTPUT FILE\n" \
  "  The program writes to standard output the adjusted feature positions" \
  " pairs, in the same format.  The coordinates {X1,Y1,X2,Y2} are usually" \
  " changed, but the radii and angles are simply copied.\n" \
  "\n" \
  "ADJUSTMENT\n" \
  "  " imft_optimize_ftpair_INFO "\n" \
  "\n" \
  "OPTIONS\n" \
  "  -image1 {FILE1} {LX1} {LY1} {HX1} {HY1}\n" \
  "  -image2 {FILE2} {LX2} {LY2} {HX2} {HY2}\n" \
  "    These mandatory arguments specify the file names and the nominal" \
  " coordinates of the low and high corners of the two images.\n" \
  "\n" \
  "  If {NX1} is the number of columns of image 1, the {X} coordinates" \
  " points on that image is mapped affinely from the range {[LX1 _ HX1]} to" \
  " (fractional) column indices in {[0 _ NX1]}, and vice-versa.\n" \
  "\n" \
  "  Similarly, if {NY1} is the number of rows of image 1, {Y} coordinates" \
  " on that image are mapped affinely from the range {[LY1 _ HY1]}" \
  " to (fractional) row indices in {[0 _ NY1]}.  And analogously for" \
  " the other image.\n" \
  "\n" \
  "  -maxAdjust {MAX_ADJ}\n" \
  "    This optional argument specifies the maximum displacement allowed" \
  " when adjusting each input point, relative tot the radius of the" \
  " comparison window ({RAD1} or {RAD2}.  For example, \"-maxAdjust 0.5\" allows" \
  " each point {(X1,Y1}) to be displaced by {RAD1/2} pixels, and each" \
  " point {(X2,Y2}) to be displaced by {RAD2/2}pixels.\n" \
  "\n" \
  "  -fix {IFIX}\n" \
  "    This optional argument specfies that the position of the feature" \
  " on image {IFIX} (1 or 2) should be kept fixed,while ajusting that of" \
  " the other. If not specified (or {IFIX} is zero), both positions are" \
  " adjusted in opposite directions.\n" \
  "\n" \
  "  -maxIter {MAX_ITER}\n" \
  "    This optional argument is the mximum number of iterations for" \
  " the non-linear optimizitation.  The default is \"-maxIter 5\".\n" \
  "\n" \
  "  -writeFeatures {OUT_PREFIX}\n" \
  "    This optional argument specifies that the features should be" \
  " extracted from the input images and written out.  Each feature is" \
  " written to a file \"{OUT_PREFIX}_ft{NNN}.png\".  The file will have" \
  " three sub-images side by side: the features extracted from images 1 and 2 at" \
  " left and right, and their difference in the middle.  Each of these" \
  " sub-images will have one pixel for each sample point used in the" \
  " comparison, multiplied by the corresponding sample weight.\n" \
  "\n" \
  "  -verbose\n" \
  "    If present, this option requests debugging printouts to {stderr}.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  Hear also.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2020-10 by J.Stolfi.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " image_match_features_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <fget.h>
#include <r2.h>
#include <r2x2.h>
#include <i2.h>
#include <ix.h>
#include <r2_aff_map.h>
#include <image_file_format.h>
#include <float_image.h>
#include <float_image_read_gen.h>
#include <float_image_aff_compare.h>
#include <float_image_aff_extract.h>
#include <float_image_write_gen.h>
#include <jsfile.h>
#include <r2_opt.h>
#include <argparser.h>

typedef struct imft_options_t 
  { char *fname[2];     /* Names of the image files. */
    r2_t L[2], H[2];    /* Nominal domain corners of each image. */
    bool_t verbose;     /* TRUE to print debugging info. */
    /* Parameters of the non-linear optimizer: */
    int32_t fix;        /* Which feature to fix (1 or 2), or 0 for none. */
    double maxAdjust;   /* Max adjustment relative to sampling disk radius. */
    int maxIter;        /* Max number of iterations. */
    /* Parameters for extracted feature output: */
    char *outPrefix;    /* The output feature file name prefix, or {NULL} if not requested. */
  } imft_options_t;
  /* Command line arguments. */
  
typedef struct imft_feature_t
  { r2_t ctr;
    double rad;
    double ang;
  } imft_feature_t;
  /* Specification of a feature on an image.  Namely the part of the image
    contained in a fuzzy disk of radius {rad} (in pixels) and center {ctr} 
    (in Client coordinates), assumed to be rotated by {ang} degrees
    counterclockwise. */
  
typedef struct imft_ftpair_t
  { imft_feature_t ft[2];
  } imft_ftpair_t;
  /* Specification of supposedly corresponding features on images 1 and 2
   (actually indexed 0 and 1).. */
    
vec_typedef(imft_ftpair_vec_t,imft_ftpair_vec,imft_ftpair_t);
 /* An extensble vector of {imft_ftpair_t}s. */

imft_options_t *imft_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

float_image_t *imft_read_image(char *fname, bool_t verbose);
  /* Reads file {fname} an image and converts its contenst to a {float_image_t}. */

imft_ftpair_vec_t imft_read_data(FILE *rd, bool_t verbose);
  /* Reads from {rd} the list of guessed matching feature pairs 
    in the format {imft_data_format_INFO}. */
  
imft_feature_t imft_read_feature(FILE *rd, bool_t verbose);
  /* Reads the specification of a single feature from a corresponding 
    feature pair as described by {imft_data_format_INFO}, namely the
    part before or after the '=' character. */

void imft_write_data(FILE *wr, imft_ftpair_vec_t *fpv, bool_t verbose);
  /* Writes to {wr} the list of adjusted matching feature pairs 
    in the format {imft_data_format_INFO}. */
 
#define imft_data_format_INFO \
  "The data file should have one feature per line, in the format\n"\
  "\n" \
  "    \"( {X1} {Y1} ) {RAD1} {ANG1} = ( {X2} {Y2} ) {RAD2} {ANG2}\"\n" \
  "\n" \
  "  The program assumes that the feature on image 1 is contained in" \
  " a fuzzy disk with approximate center {(X1,Y1)} and radius {RAD1} and" \
  " is tilted by {ANG1} degrees counterclockwise.  Likewise, it assumes" \
  " that the matching feature on image 2 is contained in a fuzzy disk" \
  " with approximate center {(X2,Y2)} and radius {RAD2}, tilted" \
  " by {ANG2} degrees.  Note that {X1,X2,Y1,Y2} are in Client" \
  " coordinates, while {RAD1,RAD1} are in pixels."

void imft_show_result(FILE *wr, r2_t ctr_ini[], r2_t ctr_adj[], r2_t pix_ini[], r2_t pix_adj[]);   
  /* Writes to {wr} the initial and final feature centers, in Client and pixel coords,
    and the total displacement in each case. */
 
void imft_optimize_ftpair
  ( imft_ftpair_t *fp,
    float_image_t *img[], 
    r2_t L[], r2_t H[],
    int32_t fix,
    double maxAdjust,
    int32_t maxIter,
    bool_t verbose
  );
  /* Adjusts the centers of the sampling disks in the two features 
    described by {fp}, on the given images. Uses {L[k],H[k]} to convert
    between Client coordinates and fractional pixel indices for image {k}. See
    {imft_optimize_ftpair_INFO}. The {verbose} flag causes diagnostics
    to be printed.*/
    
#define imft_optimize_ftpair_INFO \
  "The parameters {X1,Y1,RAD1} and {X2,Y2,RAD2} define two /sampling disks/, fuzzy" \
  " circular windows with Gaussian-like sampling weights.  The adjustment" \
  " compares the image values inside the two sampling disks, while displacing" \
  " their centers so as to minimize the discrepancy of the samples inside the disks.\n" \
  "\n" \
  "  In the current version of the program, the discrepancy is defined as the" \
  " total squared weighted differences between corresponding samples.\n" \
  "\n" \
  "  The image fragments inside the disks are implicitly scaled" \
  " by {RM/RAD1} and {RM/RAD2} (with cubic interpolation), where {RM} is" \
  " large enough for proper sampling; and rotated by {-ANG1} and {-ANG2} before" \
  " computing their discrepancy.\n" \
  "\n" \
  "  The adjustment can be applied to either {X1} or {X2}, or to both.  In any" \
  " case, the displacements are relative to the respective feature radii.  Namely," \
  " the center of the fuzzy disk on image 1 is displaced by {+dX*RAD1} pixels," \
  " and/or the center of the fuzzy disk on image 2 is displaced by {-dX*RAD2} pixels;" \
  " where {dX} is a relative adjustment in {[-MIN_ADJ _ +MIN_ADJ]}.  In Client" \
  " coordinates, the displacements to {X1} and/or {X2} are {+dX*RAD1*(HX1-LX1)/NX1} and" \
  " {-dX*RAD2*(HX2-LX2)/NX2}, respectively.  The displacements applied" \
  " to {Y1} and {Y2} have the analogous constraints, with a relative" \
  " parameter {dY} also in {[-MIN_ADJ _ +MIN_ADJ]}.\n" \
  "\n" \
  "  It follows from the above that, when both positions are being" \
  " adjusted, the displacements are opposte and equal relative to" \
  " the respective window radii. Thus, in any case, there are only" \
  " two parameters to optimize for each feature, {dX} and {dY}.  "

float_image_t *imft_read_image(char*fname, bool_t verbose);
   /* Reads an image file with given {fname} and returns the contents as a {float_image_t}.
     The image type is determined by the file name extension (".jpg", ".png", ".pnm"). */

void imft_apply_rel_adjustment
  ( r2_t ctr_ini,
    r2_t dp, 
    double rad, 
    ix_size_t sz[],
    r2_t *L, 
    r2_t *H, 
    sign_t dir, 
    r2_t *ctr_adjP,
    r2_t *pix_adjP
  );
  /* Returns in {*ctr_adjP} (if not {NULL}) the point {ctr_ini} plus the
    displacement {dir*rad*dp}. Both {ctr_ini} and {*ctr_adjP} are in
    Client coordinates, while {dir*rad*dp} is taken to be in pixels.
    Also returns in {*pix_adjP} (if not {NULL}) the point {*ctr_adjP}
    converted to fractional pixel indices.
    
    Conversion between pixel indices and Client
    coordinates is defined by the image size {sz[0..2]} (num channels,
    columns, and rows) and the Client coordinates {L,H} of the image
    corners. */

r2_t imft_pixel_from_user(r2_t *c, r2_t *L, r2_t *H, ix_size_t sz[]);
  /* Converts Client coordinates {c} of a point in to fractional pixel indices. */ 
  
void imft_write_feature_images
  ( char *prefix, 
    int32_t ip, 
    imft_ftpair_t *fp,
    r2_t L[], r2_t H[], 
    float_image_t *img[]
  );
  /* */

/* ---------------------------------------------------------------------- */
/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    fprintf(stderr, "starting...\n");
    imft_options_t *o = imft_parse_options(argc, argv);
    
    /* Read the images: */
    fprintf(stderr, "reading images...\n");
    float_image_t *img[2];
    img[0] = imft_read_image(o->fname[0], o->verbose);
    img[1] = imft_read_image(o->fname[1], o->verbose);
    
    /* Read the point pairs: */
    fprintf(stderr, "reading feature specs...\n");
    imft_ftpair_vec_t fpv = imft_read_data(stdin, o->verbose);
    int np = fpv.ne; /* Number of point pairs. */
    
    for (uint32_t ip = 0;  ip < np; ip++)
      { fprintf(stderr, "optimizing feature %d...\n", ip);
        imft_ftpair_t *fp = &(fpv.e[ip]);
        imft_optimize_ftpair
          ( fp,
            img, o->L, o->H,
            o->fix,
            o->maxAdjust,
            o->maxIter,
            o->verbose
          );

        if (o->outPrefix != NULL)
          { fprintf(stderr, "writing the feature pair image...\n");
            imft_write_feature_images(o->outPrefix, ip, fp, o->L, o->H, img);
          }
      }
      
    fprintf(stderr, "writing adjusted specs...\n");
    imft_write_data(stdout, &fpv, o->verbose);
    fflush(stdout);
    return 0;
  }

void imft_optimize_ftpair
  ( imft_ftpair_t *fp,
    float_image_t *img[], 
    r2_t L[], r2_t H[],
    int32_t fix,
    double maxAdjust,
    int32_t maxIter,
    bool_t verbose
  )
  {
    bool_t debug = FALSE;
    
    auto double sqr_mismatch(int nit, r2_t dpt[], i2_t isct);
      /* Computes the mismatch for relative displacement {dpt[0]}. 
        Requires {nit} to be 1.  For now, {isct} should be {(0,0)}. */
  
    r2_t ctr_ini[2];  /* Original centers (Client coords): */
    r2_t pix_ini[2];  /* Original centers (Pixel coords): */
    r2_aff_map_t aff[2];/* Affine maps to use for comparison: */
    double dtor = M_PI/180; /* Degree to radian conversion. */
    for (uint32_t k = 0;  k < 2; k++)
      { ctr_ini[k] = fp->ft[k].ctr;
        pix_ini[k] = imft_pixel_from_user(&(ctr_ini[k]), &(L[k]), &(H[k]), img[k]->sz);
        aff[k] = r2_aff_map_rot_scale(fp->ft[k].ang*dtor, fp->ft[k].rad);
      }
        
    /* Optimize {dp[0]}: */
    i2_t iscale = (i2_t){{ 0,0 }}; /* For now. */
    r2_t arad[1]; arad[0] = (r2_t){{ maxAdjust, maxAdjust }}; /* Max rel displacement. */ 
    r2_t astp[1]; astp[0] = (r2_t){{ maxAdjust/10, maxAdjust/10 }}; /* Req accuracy. */
    r2_t dp[1]; /* The relative displacement for the sampling windows. */
    double sqd; /* The square mismatch for {dp}. */
    r2_opt_single_scale_quadopt
      ( 1, iscale,
        &sqr_mismatch,
        arad, astp,
        dp, 
        &(sqd),
        debug
      );
    
    /* Update the feature pair with the optimal displacement {dp[0]}: */
    r2_t ctr_adj[2];
    r2_t pix_adj[2];
    for (uint32_t k = 0;  k < 2; k++)
      { sign_t dir = (sign_t)(1 - 2*k);
        if (fix != k+1)
          { imft_apply_rel_adjustment
              ( ctr_ini[k], dp[0], fp->ft[k].rad, 
                img[k]->sz, &(L[k]), &(H[k]), 
                dir, 
                &(ctr_adj[k]), &(pix_adj[k])
              );
            fp->ft[k].ctr = ctr_adj[k];
          }
        else
          { fp->ft[k].ctr = ctr_ini[k]; }
      }
      
    fprintf(stderr, "reporting the result...\n");
    imft_show_result(stderr, ctr_ini, ctr_adj, pix_ini, pix_adj);

    return;
    
    /* Internal implementations: */

    double sqr_mismatch(int nit, r2_t dpt[], i2_t isct)
      { assert(nit == 1);
        assert((isct.c[0] == 0) && (isct.c[1] == 0)); /* For now. */
        /* Apply {dpt[0]} to the original centers to obtain the tentative ones: */
        for (uint32_t k = 0;  k < 2; k++)
          { r2_t pixt;
            sign_t dir = (sign_t)(1 - 2*k);
            if (fix != k+1)
              { imft_apply_rel_adjustment
                  ( ctr_ini[k], dpt[0], fp->ft[k].rad, img[k]->sz, &(L[k]), &(H[k]), dir, NULL, &(pixt) );
                /* Insert the new centers in the affine maps: */
              }
            else
              { pixt = pix_ini[k]; }
            aff[k].disp = pixt;
          }
        /* Evaluate the discrepancy: */
        double d2 = float_image_aff_compare(img[0], &(aff[0]), img[1], &(aff[1]), NULL, NULL);
        return d2;
      }
  }

imft_ftpair_vec_t imft_read_data(FILE *rd, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- enter imft_read_data ---\n"); }
    imft_ftpair_vec_t fpv = imft_ftpair_vec_new(20); /* To be expanded/trimmed. */
    int np = 0; /* Number of data pairs. */

    while (TRUE)
      { /* Discards blank or comment lines: */
        if (fget_test_comment_or_eol(rd, '#', NULL)) { continue; }
        if (feof(rd)) { break; }
        imft_ftpair_vec_expand(&fpv, np);
        imft_ftpair_t *fp = &(fpv.e[np]);
        fp->ft[0] = imft_read_feature(rd, verbose);
        fget_skip_spaces_and_match(rd, "=");
        fp->ft[1] = imft_read_feature(rd, verbose);
        fget_comment_or_eol(rd, '#', NULL);
        np++;
      }
    imft_ftpair_vec_trim(&fpv, np);
    if (verbose) { fprintf(stderr, "--- exit imft_read_data ---\n"); }
    return fpv;
  }

imft_feature_t imft_read_feature(FILE *rd, bool_t verbose)
  { 
    imft_feature_t ft;
    fget_skip_spaces_and_match(rd, "(");
    ft.ctr.c[0] = fget_double(rd);
    ft.ctr.c[1] = fget_double(rd);
    fget_skip_spaces_and_match(rd, ")");
    ft.rad = fget_double(rd);
    ft.ang = fget_double(rd);
    return ft;
  }

void imft_write_data(FILE *wr, imft_ftpair_vec_t *fpv, bool_t verbose)
  {
    int32_t np = fpv->ne;
    for (uint32_t ip = 0;  ip < np; ip++)
      { imft_ftpair_t *fp = &(fpv->e[ip]);
        for (uint32_t k = 0;  k < 2; k++)
          { imft_feature_t *ft = &(fp->ft[k]);
            fprintf(wr, "( %.6f %.6f )", ft->ctr.c[0], ft->ctr.c[1]);
            fprintf(wr, " %.3f %.3f", ft->rad, ft->ang);
            if (k == 0) { fprintf(wr, " = "); }
          }
        fprintf(wr, "\n");
      }
    fflush(wr);
  }

void imft_show_result(FILE *wr, r2_t ctr_ini[], r2_t ctr_adj[], r2_t pix_ini[], r2_t pix_adj[])
  {
    for (uint32_t k = 0;  k < 2; k++)
      { fprintf(wr, "adjusted center on image %d:\n", k+1);
        for (uint32_t sys = 0;  sys < 2; sys++) 
          { fprintf(wr, "  %s:", (sys == 0 ? "Client" : "Pixels"));
            r2_t ini = (sys == 0 ? ctr_ini[k] : pix_ini[k]);
            r2_t adj = (sys == 0 ? ctr_adj[k] : pix_adj[k]);
            r2_t dif; r2_sub(&adj, &ini, &dif);
            fprintf(wr, " ( %12.6f %12.6f )", ini.c[0], ini.c[1]);
            fprintf(wr, " + ( %12.6f %12.6f )", dif.c[0], dif.c[1]);
            fprintf(wr, " = ( %12.6f %12.6f )", adj.c[0], adj.c[1]);
            fprintf(wr, "\n");
          }
      }
    fprintf(wr, "\n");
    fflush(wr);
  }

float_image_t *imft_read_image(char*fname, bool_t verbose)
  { image_file_format_t ffmt = image_file_format_from_name(fname);
    float_image_t *img = float_image_read_gen_named
      ( fname, ffmt, 0.0f, 1.0f, NULL, NULL, NULL, verbose);
    return img;
  }
    
void imft_apply_rel_adjustment
  ( r2_t ctr_ini,
    r2_t dp, 
    double rad, 
    ix_size_t sz[],
    r2_t *L, 
    r2_t *H, 
    sign_t dir, 
    r2_t *ctr_adjP, 
    r2_t *pix_adjP
  )
  {
    r2_t ctr_adj;
    for (uint32_t j = 0;   j < 2; j++)
      { /* Scale factor to convert pixel indices to Client coordinates: */
        double ptoc = (H->c[j] - L->c[j])/((double)sz[j+1]);
        /* Grab original Client coordinate: */
        /* Compute absolute displacement in pixels: */
        double d_pix = dir*dp.c[j]*rad;
        /* Convert to displacement in Client coordinates: */
        double d_ctr = d_pix*ptoc;
        /* Compute the Client coords: */
        ctr_adj.c[j] = ctr_ini.c[j] + d_ctr;
      }
    if (ctr_adjP != NULL) { (*ctr_adjP) = ctr_adj; }
    if (pix_adjP != NULL) { (*pix_adjP) = imft_pixel_from_user(&(ctr_adj), L, H, sz); }

  }
  
r2_t imft_pixel_from_user(r2_t *ctr, r2_t *L, r2_t *H, ix_size_t sz[])
  {
    r2_t pix;
    for (uint32_t j = 0;   j < 2; j++)
      { /* Scale factor to convert pixel indices to Client coordinates: */
        double ptoc = (H->c[j] - L->c[j])/((double)sz[j+1]);
        /* Grab original Client coordinate: */
        pix.c[j] = (ctr->c[j] - L->c[j])/ptoc;
      }
    return pix;
  }

void imft_write_feature_images
  ( char *prefix, 
    int32_t ip, 
    imft_ftpair_t *fp,
    r2_t L[], r2_t H[], 
    float_image_t *img[]
  )
  {
    int32_t NC = (int32_t)img[0]->sz[0];
    demand(NC == (int32_t)img[1]->sz[0], "inconsistent channel counts");
    
    /* Construct the affine maps: */
    double dtor = M_PI/180.0; /* Convert degrees to radians. */
    r2_aff_map_t aff[2];
    for (uint32_t k = 0;  k < 2; k++) 
      { r2_t ctrk = fp->ft[k].ctr;
        r2_t pixk = imft_pixel_from_user(&(ctrk), &(L[k]), &(H[k]), img[k]->sz);
        aff[k] = r2_aff_map_rot_scale(fp->ft[k].ang*dtor, fp->ft[k].rad);
        aff[k].disp = pixk;
      }
    
    /* Extract the features: */
    r2_t step;
    i2_t size;
    (void)float_image_aff_compare(img[0], &(aff[0]), img[1], &(aff[1]), &step, &size); 
    fprintf(stderr, "extracted feature size = (%d %d)\n", size.c[0], size.c[1]);
    float_image_t *ftr[2];
    int32_t NXF = size.c[0];
    int32_t NYF = size.c[1];
    float_image_t *res = float_image_new(NC, 3*NXF, NYF);
    for (uint32_t k = 0;  k < 2; k++)
      { ftr[k] = float_image_aff_extract(img[k], &(aff[k]), step, size);
        int32_t xminA = 2*NXF*k, xmaxA = xminA + NXF - 1;
        float_image_assign_rectangle(res, xminA, xmaxA, 0,NYF-1, ftr[k], 0, 0);
      }
      
    /* Compute the difference image: */
    for (uint32_t iy = 0;  iy < NYF; iy++)
      { for (uint32_t ix = 0;  ix < NXF; ix++)
          { for (uint32_t ic = 0;  ic < NC; ic++)
              { float v0 = float_image_get_sample(ftr[0], ic, ix, iy);
                float v1 = float_image_get_sample(ftr[1], ic, ix, iy);
                float vdif = (float)(0.5 + (v1 - v0)/2);
                float_image_set_sample(res, ic, ix + NXF, iy, vdif);
              }
          }
      }
    for (uint32_t k = 0;  k < 2; k++) { float_image_free(ftr[k]); ftr[k] = NULL; }
    
    /* Write the composite image: */
    char *fname = jsprintf("out/%s_%03d.png", prefix, ip);
    image_file_format_t ffmt = image_file_format_PNG;
    fprintf(stderr, "writing to file %s ...\n", fname);
    float_image_write_gen_named(fname, res, ffmt, 0.0, 1.0, NAN, NAN, FALSE);
    free(fname);
    float_image_free(res); res = NULL;
    
  }

imft_options_t *imft_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    imft_options_t *o = (imft_options_t *)notnull(malloc(sizeof(imft_options_t)), "no mem"); 

    /* Parse keyword-based arguments: */
    for (uint32_t k = 0;  k < 2; k++)
      { char *key = (k == 0 ? "-image1" : "-image2" );
        argparser_get_keyword(pp, key);
        o->fname[k] = argparser_get_next(pp);
        o->L[k].c[0] = argparser_get_next_double(pp, -100000, +100000);
        o->L[k].c[1] = argparser_get_next_double(pp, -100000, +100000);
        o->H[k].c[0] = argparser_get_next_double(pp, -100000, +100000);
        o->H[k].c[1] = argparser_get_next_double(pp, -100000, +100000);
      }
    
    if (argparser_keyword_present(pp, "-fix"))
      { o->fix = (int32_t)argparser_get_next_int(pp, 1, 2); }
    else
      { o->fix = 0; }

    if (argparser_keyword_present(pp, "-maxIter"))
      { o->maxIter = (int32_t)argparser_get_next_int(pp, 0, 100000); }
    else
      { o->maxIter = 5; }

    if (argparser_keyword_present(pp, "-maxAdjust"))
      { o->maxAdjust = argparser_get_next_double(pp, 0, 100.0); }
    else
      { o->maxAdjust = 0.0; }

    if (argparser_keyword_present(pp, "-writeFeatures"))
      { o->outPrefix = argparser_get_next(pp); }
    else
      { o->outPrefix = NULL; }

    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

vec_typeimpl(imft_ftpair_vec_t,imft_ftpair_vec,imft_ftpair_t); 
