#define PROG_NAME "image_stitch_n"
#define PROG_DESC "Finds projective map between N images, given corresponding points"
#define PROG_VERS "1.0"

// Last edited on 2023-10-09 19:36:19 by stolfi

#define image_stitch_n_C_COPYRIGHT \
    "© 2002 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    { -image {ID} {LX} {LY} {HX} {HY} }.. \\\n" \
  "    { -match {IDA} {IDB} {PTPAIRFILE} }.. \\\n" \
  "    -outPrefix {STRING} \\\n" \
  "    [ -perpective {PERSPBOOL} ] \\\n" \
  "    [ -maxIter {MAXITER} ] \\\n" \
  "    [ -maxErr {MAXERR} ] \\\n" \
  "    [ -eqDistortion {EQDISTBOOL} ] \\\n" \
  "    [ -minScale {MINSCALE} ] [ -maxScale {MAXSCALE} ] \\\n" \
  "    [ -eqRotation {EQROTBOOL} ] \\\n" \
  "    [ -center {XCTR} {YCTR} | -minCoords {XMIN} {YMIN} ] \\\n" \
  "    [ -verbose {VBOOL} ] \\\n" \
  "    > {MATRIXFILE} \\\n"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "   Reads lists of pairs of corresponding points in two or more images.  Computes" \
  " an affine or projective map for each image that produce a good match for" \
  " those point pairs.  Optionally, it may optimize the overall" \
  " distortion, scaling, rotation and position of the mapped images.\n" \
  "\n" \
  "   Each image may have its own coordinate system, with arbitrary" \
  " length units.  This coordinate system is used to specify the coordinates" \
  " of points in each image, the image's nominal domain, and its map's matrix.\n" \
  "\n" \
  "  Each image has a /domain/ which is a user-specified convex" \
  " quadrilateral.  The intersection of the two diagonals is called" \
  " the /centroid/ of the domain.  Note that the centroid of a projectively" \
  " mapped domain is the projective image of the original centroid.\n" \
  "\n" \
  "   The projective maps used by this program are constrained" \
  " to map each image's domain to a finite (and not very" \
  " large) quadrilateral.  Therefore this  program" \
  " is not suitable for stitching sets of panoramic photos" \
  " whose combined span exceeds (or is even close to) 180 degrees.\n" \
  "\n" \
  "INPUT FILES\n" \
  "   Each point pair file (specified by the \"-match\" option) relates" \
  " two images.  It must have one point pair per" \
  " line, in the format \"( {X1} {Y1} ) = ( {X2} {Y2} )\" where" \
  " {(X1,Y1)} is a point in the domain of one image, and {(X2,Y2)} is" \
  " a point in the domain of the other image.  Blank lines and '#'-comments are ignored.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "   Each projective map is written out as a pair of" \
  " files \"{PREFIX}-matrix-{ID}-{MDIR}.txt\" where {PREFIX} is" \
  " specified by the \"-outPrefix\" argument, {ID} is the image label specified" \
  " by the \"-image\" argument, and {MDIR} is either \"dir\" or \"inv\" for the" \
  " direct and inverse homogeneous matrix, respectively.  Each matrix is written as nine" \
  " numbers in row-by-row order. \n" \
  "\n" \
  "   For each point pair list specified by a \"-match\" argument, the" \
  " program also writes a file called \"{PREFIX}-match-{IDA}-{IDB}.txt\" that" \
  " contains the given point pairs mapped by the respective image maps, one per" \
  " line, as four blank-separated numbers. \n" \
  "\n" \
  "OPTIONS\n" \
  "  For the boolean options below," \
  " \"t\", \"T\", \"true\", \"TRUE\", \"y\", \"Y\", \"yes\", \"YES\", or" \
  " 1 mean true, and (\"f\", \"F\", \"false\", \"FALSE\"," \
  " \"n\", \"N\",\"no\", \"NO\", or" \
  " 0 mean false.\n" \
  "\n" \
  "  -image {ID} {X0} {Y0}  {X1} {Y1}  {X2} {Y2}  {X3} {Y3}\n" \
  "    This flag must be given once for each image.  It" \
  " specifies an arbitrary label {ID} and the nominal" \
  " domain for the image, namely the convex quadrilateral with corners" \
  " {P0=(X0,Y0)}, {P1=(X1,Y1)}, {P2=(X2,Y2)} and {P3=(X3,Y3)}" \
  " which must be given in counterclockwise order relatve to the" \
  " image's coordinate system. Typically these are" \
  " {P0=(0,0)}, {P1=(W,0)}, {P2=(W,H)} and {P3=(0,H)}," \
  " where {W} is the width of image in pixels, and" \
  " {H} is its height in pixels.\n" \
  "\n" \
  "    The order of the images is important: every image must have at least 3 points" \
  " matched to points in the previous images in the list.  Moreover, unless explicitly" \
  " specified otherwise, the first image will remain undeformed while" \
  " the remaining ones will be deformed to match it.\n" \
  "\n" \
  "    The image domain is used for distortion, scale, rotation and" \
  " centering purposes.  It is also used to validate the points" \
  " in the matched point pair files, which must all be" \
  " inside this quadrilateral. (Because of arithmetic roundoff" \
  " errors, this test may fail if the data points are inside the domain but too close" \
  " to its sides. In that case one may have to fudge the domain corners and/or the data points" \
  " to remove this bogus error.)\n" \
  "\n" \
  "  -match {IDA} {IDB} {PTPAIRFILE}\n" \
  "    This flag may be given (at most once) for each pair of distinct images," \
  " identified by the labels {IDA} and {IDB} (which must be distinct).  The file {PTPAIRFILE} must" \
  " then contain one or more point pairs, in the format specified in the INPUT FILES section.\n" \
  "\n" \
  "  The point pair files pairs must be enough to connect all images.  There may be at" \
  " most one \"-match\" argument for any pair of images {IDA,IDB}.\n" \
  "\n" \
  "  -outPrefix {STRING}\n" \
  "    This mandatory argument is the prefix for the file names of all output files.\n" \
  "\n" \
  "  -perpective {PERSPBOOL}\n" \
  "    If the boolean value {PERSPBOOL} is false, the" \
  " program will use only affine maps.  If {PERSPBOOL} is true," \
  " it may also use perspective distortion maps.  The default is \"-perspective YES\".\n" \
  "\n" \
  "  -maxIter {MAXITER}\n" \
  "    This optional argument is the number of iterations for" \
  " the non-linear optimization of individual maps (including individual" \
  " perspective distortion, if allowed).  If {MAXITER} is zero, the non-linear optimization" \
  " is skipped, and all maps will be affine (except for a global perspective" \
  " map, if \"-eqDistortion\" is requested).  The default is \"-maxIter 5\".\n" \
  "\n" \
  "  -maxErr {MAXERR}\n" \
  "    This optional argument is the convergence criterion for non-linear" \
  " optimization.  The optimization will stop when the root mean square" \
  " mismatch of the mapped points (after all distorion and" \
  " scale equalizations) is less than {MAXERR}.  If {MAXERR} is zero" \
  " (the default), the optimization will continue for {MAXITER} iterations.\n" \
  "\n" \
  "  -eqDistortion {EQDISTBOOL}\n" \
  "    If the boolean value {EQDISTBOOL} is true, the" \
  " image maps will be adjusted by a single non-uniform scaling" \
  " and perspective map so that the domain" \
  " shape distortion is equally distributed among all" \
  " images.  If {EQDISTBOOL} is false," \
  " the domain of the first image will remain" \
  " undistorted (except perhaps for rotation, translation" \
  " and scale), and the rest will be distorted as needed to" \
  " match it.  The default is \"-eqDistortion NO\".\n" \
  "\n" \
  "  -minScale {MINSCALE}\n" \
  "  -maxScale {MAXSCALE}\n" \
  "    If any of these arguments is specified, nonzero, and finite," \
  " image maps will be adjusted by a global uniform scaling (after distortion" \
  " equalization, if any) so that the linear scaling factor at the corners of every" \
  " domain is at most {MAXSCALE} and/or at least {MINSCALE}.  If both are" \
  " specified but impossible to satisfy at the same time, the program will issue" \
  " a warning and ignore them.\n" \
  "\n" \
  "  -eqRotation {EQROTBOOL}\n" \
  "    If the boolean value {EQROTBOOL} is true, the" \
  " image maps will be adjusted by a single global rotation map (after distortion" \
  " and scale equalization, if requested) so that the mean rotation" \
  " angle is evenly distributed among all images.  If {EQROTBOOL} is false, the" \
  " domain of the first image will remain unrotated, and the rest" \
  " will be rotated as needed to match it.  The default is \"-eqRotation NO\".\n" \
  "\n" \
  "  -center {XCTR} {YCTR}\n" \
  "  -minCoords {XMIN} {YMIN}\n" \
  "    These arguments specify a global translation to be applied to all" \
  " image maps (after distortion," \
  " scale and rotation equalization, if any).  If \"-center\" is specified" \
  " and {XCTR,YCTR} are finite, the" \
  " translation will be such that the bounding box of all mapped domains will be centered" \
  " at {(XCTR,YCTR)}. If \"-minCoords\" is specified and {XMIN,YMIN} is" \
  " finite, the translation will be such that the minimum {X} coordinate" \
  " of the mapped domains is {XMIN}, and the minimum {Y} coordinate" \
  " is {YMIN}.  These two options are mutually exclusive.  If neither" \
  " option is requested, the" \
  " domain of the first image will remain at its original" \
  " position, and the rest" \
  " will be translated as needed to match it.\n" \
  "\n" \
  "  -verbose {VBOOL}\n" \
  "    If {VBOOL} is true, progress and debugging information" \
  " will be issued to {stderr}.  The default is \"-verbose NO\".\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  Hear also.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 01/jan/2012 by J.Stolfi.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " image_stitch_n_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <hr2.h>
#include <r2.h>
#include <r2x2.h>
#include <r3.h>
#include <r3x3.h>
#include <rn.h>
#include <jsfile.h>
#include <fget.h>
#include <interval.h>
#include <interval_io.h>
#include <sve_minn.h>
#include <argparser.h>

/* !!! Make sure that there is at most one -match for each pair of images !!! */

/* ---------------------------------------------------------------------- */
/* COMMAND LINE OPTIONS */

typedef struct options_image_t /* Command line info about one image. */
  { char *id;           /* Image label. */
    r2_t P[4];          /* Nominal domain corners of image. */
  } options_image_t;
vec_typedef(options_image_vec_t,options_image_vec,options_image_t);

typedef struct options_match_t  /* Command line info about one point pair list. */
  { char *idA;          /* Image A label. */
    char *idB;          /* Image B label. */
    char *fname;        /* Name of point pairs file. */
  } options_match_t;
vec_typedef(options_match_vec_t,options_match_vec,options_match_t);

typedef struct options_t  /* Command line arguments. */
  { options_image_vec_t image;    /* Image-related options. */
    options_match_vec_t match;    /* Corresponding point pair options. */
    bool_t verbose;               /* TRUE to print debugging info. */
    char *outPrefix;              /* Prefix for output filenames. */
    /* Mapping and equalization parameters: */
    bool_t perspective;      /* Allow perspective distortion. */
    bool_t eqDistortion;     /* Perform distortion equalization. */
    double minScale;         /* Min linear scale at domain corners. */
    double maxScale;         /* Max linear scale at domain corners. */
    bool_t eqRotation;       /* Perform rotation equalization. */
    r2_t center;             /* Center of the mapped domains. */
    r2_t minCoords;          /* Minimum coordinates of the mapped domains. */
    /* Parameters of the non-linear optimizer: */
    int maxIter;                  /* Max number of iterations. */
    double maxErr;                /* Convergence criterion. */
  } options_t;

options_t *imsn_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

/* ---------------------------------------------------------------------- */
/* MAIN PROCEDURES */

typedef struct match_data_t /* A list of matched point pairs. */
  { int ixA;        /* Image A index. */
    r2_vec_t pA;    /* Points in image A domain. */
    int ixB;        /* Image B index. */
    r2_vec_t pB;    /* Points in image B domain. */
  } match_data_t;
vec_typedef(match_data_vec_t,match_data_vec,match_data_t);

int main(int argc, char **argv);

void imsn_read_all_data
  ( options_image_vec_t *image, /* (IN) Image options. */
    options_match_vec_t *match, /* (IN) Point pair options. */
    bool_t verbose,             /* (IN) Print diagnostics. */
    match_data_vec_t *mt        /* (OUT) Point pair lists. */  
  );
  /* Reads the lists of point pairs {mt.e[0..nmt-1]},
    specified in {match.e[0..nmt-1]}, where {nmt = match.ne}.
    The vector {mt} is (re)allocated and/or trimmed by the
    procedure to have size {nmt}. */

vec_typedef(hr2_pmap_vec_t,hr2_pmap_vec,hr2_pmap_t);

void imsn_compute_initial_affine_pmaps
  ( options_image_vec_t *image,   /* (IN) Image options. */      
    match_data_vec_t *mt,         /* (IN) Point pair options. */ 
    bool_t verbose,               /* (IN) Print diagnostics. */  
    hr2_pmap_vec_t *M             /* (OUT) Affine maps. */  
  );
  /* Computes affine maps {M.e[0..nim-1]} that map each image to a
    common mosaic, so that the match information {mt.e[0..mt.ne-1]} is
    best satisfied; where {nim = image.ne}. 
    
    The vector {M} is (re)allocated and/or trimmed by the procedure to
    have size {nim}. The computed maps will be affine; that is, their
    elements {[0][0]}, {[1][0]} and {[2][0]} will be 1, 0, 0. The
    first map {M.e[0]} will be the identity. */

void imsn_equalize_pmap_distortion
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M,
    bool_t perspective,
    bool_t verbose
  );
  /* Applies to all pmaps {M.e[0..nim-1]} a single distortion
    map that minimizes the total distortion of all image domains.
    If {perspective} is true, uses a general perspective distortion map
    with fixed uniform scaling, rotation and translation component.
    If {perspective} is false, uses an affine map with 
    only shear and non-uniform scaling. */

void imsn_adjust_global_scale
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M,
    double minScale,
    double maxScale,
    bool_t verbose
  );
  /* Applies to all pmaps {M.e[0..nim-1]} a single uniform scaling map
    so as to satisfy the constraints specified by {minScale,maxScale}.
    
    If {minScale} (resp. {maxScale}) is not zero, the linear
    magnification at the domain corners will be adjusted if needed so
    that it is at least {minScale} (resp. at most {maxScale}). If
    {minScale} and {maxScale} are both nonzero but impossible to honor
    at the same time, they are ignored. */

void imsn_optimize_all_pmaps
  ( options_image_vec_t *image,
    match_data_vec_t *mt,
    hr2_pmap_vec_t *M, 
    bool_t perspective,
    int maxIter,
    double maxErr, 
    bool_t verbose
  );
  /* Optimizes the projective maps {M[0..nim-1]} that map each 
    image to a common mosaic so that the match information
    {mt[0..nmt-1]} is best satisfied, the deformation 
    of the image domains is spread out somewhat evenly
    among the images, and no image has its resolution reduced.
    Here {nim = image.ne} and {nmt = mt.ne}.  */

void imsn_equalize_rotation
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M, 
    bool_t verbose
  );
  /* Applies to all pmaps {M.e[0..nim-1]} a single rotation
    map that equalizes the mean rotation angle at the centroid
    of each domain. */

void imsn_adjust_translation
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M,
    r2_t *center,
    r2_t *minCoords,
    bool_t verbose
  );
  /* Applies to all pmaps {M.e[0..nim-1]} a single translation.
    If {center} is finite, the cener of the boundin box of all mapped domains 
    will end up at {center}.  If {minCoords} is finite,
    the min {X} and min {Y} coordinates of the mapped domains
    will be {minCoords}. */

void imsn_show_all_pmaps
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M, 
    char *descr
  );
  /* Writes (to stderr) the direct and inverse matrices of {M.e[0..nim-1]}, 
    each labeled by the corresponding matrix {id} and by the title {descr}. */

void imsn_output_all_pmaps
  ( char *outPrefix,
    options_image_vec_t *image,
    hr2_pmap_vec_t *M, 
    bool_t verbose
  );
  /* Writes the direct and inverse matrices of {M.e[0..nim-1]} to 
    text files called "{outPrefix}-{id}-dir.txt" and "{outPrefix}-{id}-inv.txt"
    where {id} is the matrix {id}. */

void imsn_show_all_matches
  ( options_image_vec_t *image,
    match_data_vec_t *mt, 
    hr2_pmap_vec_t *M, 
    char *descr
  );
  /* Writes (to stderr) the points in each point pair list {mt.e[0..mt.ne-1]}
    mapped by their respective projective maps in {{M.e[0..nim-1]},
    showing the point differences and their rms and maximum values;
    all labeled by the corresponding matrix {id}s and by the title {descr}. */

void imsn_output_all_matches
  ( char *outPrefix,
    options_image_vec_t *image,
    match_data_vec_t *mt, 
    hr2_pmap_vec_t *M, 
    bool_t verbose
  );
  /* Writes each point pair list {mt.e[0..mt.ne-1]}
    mapped by their respective projective maps in {{M.e[0..nim-1]}
    to files "{outPrefix}-match-{idA}-{idB}.txt", where {idA} and {idB}
    are the corresponding image identifiers. */

/* ---------------------------------------------------------------------- */
/* LOWER-LEVEL PROCEDURES */

void imsn_read_point_pair_list(char *fname, r2_vec_t *pA, r2_vec_t *pB, bool_t verbose);
/* Reads from file {fname} a list of point pairs {pA[i],pB[i]},
  where {pA[i]} is in image A and {pB[i]} is in image B.
  Each pair should be in the format "( hA vA ) = ( hB vB )" where 
  {hA,hB} are column indices and {vA,vB} are row indices,
  both counted from 0.  The vectors {pA,pB} are (re)allocated
  and/or trimmed by the procedure as needed. */

int imsn_lookup_image_id(char *id, options_image_vec_t *image);
  /* Returns the index of {k} the entry {image.e[k]} with the given {id}
    field. Bombs out if there is no such entry. */

bool_t imsn_check_point_containment(options_image_t *im, r2_vec_t *p);
  /* Returns true iff all points in the list {p} are inside the 
    domain of image {im}.   */

void imsn_collect_previous_matches
  ( int m, 
    hr2_pmap_vec_t *M, 
    match_data_vec_t *mt, 
    r2_vec_t *pA, 
    r2_vec_t *pB,
    bool_t verbose
  );
  /* Returns in {*pA} the list of all points of images {0..m-1} that 
    are matched by {mt} with points of image {m}, and in {*pB}
    the corresponding points of image {m}.  The points in {pA}
    (but not those in {pB}) are mapped by the corresponding
    matrices from {M.e[0..m-1]}. The vectors {*pA,*pB}
    are (re) allocated and trimmed as needed. */

void imsn_all_mapped_domains_bbox(options_image_vec_t *image, hr2_pmap_vec_t *M, interval_t B[]);
  /* Returns in {B[0]} and {B[1]} the range of X and Y coordinates of the domains of the images
    {image.e[0..nim-1]} mapped by {M.e[0..nim-1]}. */
    
void imsn_mapped_domain_bbox(options_image_t *im, hr2_pmap_t *M, interval_t B[]);
  /* Returns in {B[0]} and {B[1]} the range of X and Y coordinates of the domain of image
    {im} mapped by {M}. */

void imsn_translate_all_images(hr2_pmap_vec_t *M, r2_t *vec);
  /* Compose all image maps {M.e[0..nim-1]} with a translation by {vec}. */

hr2_pmap_t imsn_compute_affine_pmap(r2_vec_t *pA, r2_vec_t *pB, bool_t verbose);
  /* Computes an affine map {M} such that {pB} mapped by {M} 
    matches {pA}, approximately. */

void imsn_map_point(r2_t *p, r3x3_t *M, r2_t *q);
  /* Maps the Cartesian point {p} through the projective map with 
    homogeneous matrix {M}, stores result in {*q}. Returns {(INF,INF)} 
    if the resulting point is at infinity or beyond. */

void imsn_bar(r2_vec_t *p, r2_t *bar);
  /* Computes the barycenter {bar} of the point set {*p}. */

double imsn_mean_err_sqr(r2_vec_t *pA, r3x3_t *MA, r2_vec_t *pB, r3x3_t *MB);
  /* Cmputes the mean squared error between the points of {pA} mapped by 
    {MA} and the corresponding points of {pB} mapped by {MB}. */

double imsn_deform_sqr(r2_t P[], r3x3_t *M);
  /* Computes a quadratic deformation term for the quadrilateral with  
    corners {P[0..3]} when mapped through the projective map {M}. 
    The result is 0 for Euclidean maps, and positive for other maps. */

r2_t imsn_domain_centroid(r2_t Q[]);
  /* Returns the centroid of the convex quadrilateral with corners
    { Q[0..3]}, namely the intersection of the two diagonals. */

void imsn_show_pmap(hr2_pmap_t *M, char *id, char *descr);
  /* Writes to stderr the direct and inverse matrices of {M}, labeled 
    by the matrix id {id} and by the title {descr}. */

void imsn_show_mapped_domain(options_image_t *im, hr2_pmap_t *M);
  /* Writes to stderr the domain of image {im} mapped by {M}, as well as its centroid. */

void imsn_show_bbox_and_center(interval_t B[], char *what);
  /* Writes to stderr the bounding box {B[0..1]}, as well as its center.
    The text {what} is preixed to the line. */

void imsn_output_pmap(char *outPrefix, char *id, hr2_pmap_t *M, bool_t verbose);
  /* Writes to files called "{outPrefix}-matrix-{id}-dir.txt" and "{outPrefix}-matrix-{id}-inv.txt"
    the direct and inverse matrices of {M}. */

void imsn_output_matrix(char *outPrefix, char *id, char *dir, r3x3_t *M, bool_t verbose);
  /* Writes the matrix {M} to a file called "{outPrefix}-matrix-{id}-{dir}.txt". */

void imsn_show_match
  ( r2_vec_t *pA, 
    r3x3_t *MA, 
    char *idA, 
    r2_vec_t *pB, 
    r3x3_t *MB, 
    char *idB, 
    char *descr
  );
  /* Writes (to stderr) the points of {pA} from image {idA} mapped by
    {MA}, compared to the corresponding points of {pB} from image
    {idB} mapped by {MB}, showing the point differences and their rms
    and maximum values, all with the title {descr}. */

void imsn_output_match
  ( char *outPrefix,
    r2_vec_t *pA, 
    r3x3_t *MA, 
    char *idA, 
    r2_vec_t *pB, 
    r3x3_t *MB, 
    char *idB,
    bool_t verbose
  );
  /* Writes to file called "{outPrefix}-match-{idA}-{idB}.txt"
    the points {pA} maped by {MA] and the points {pB} mapped by {MB}. */

/* ---------------------------------------------------------------------- */
/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    options_t *o = imsn_parse_options(argc, argv);
    
    int nim = o->image.ne; /* Number of images. */
    int nmt = o->match.ne; /* Number of point pair lists. */
    
    /* Read all matched point pair lists {mt.e[0..nmt-1]}. */
    match_data_vec_t mt = match_data_vec_new(0);
    imsn_read_all_data(&(o->image), &(o->match), o->verbose, &mt);
    assert(mt.ne == nmt);
    
    /* Compute the initial affine maps {M0.e[0...ne]}. */
    hr2_pmap_vec_t M = hr2_pmap_vec_new(0);
    imsn_compute_initial_affine_pmaps(&(o->image), &mt, o->verbose, &M);
    assert(M.ne == nim);
                                        
    if (o->verbose)
      { imsn_show_all_pmaps(&(o->image), &M, "affine");
        imsn_show_all_matches(&(o->image), &mt, &M, "affine");
      }

    /* Equalize the distortion among all image maps, if requested: */
    if (o->eqDistortion)
      { imsn_equalize_pmap_distortion(&(o->image), &M, o->perspective, o->verbose); }
        
    /* Adjust the magnification of all image maps, as  requested: */
    imsn_adjust_global_scale(&(o->image), &M, o->minScale, o->maxScale, o->verbose);
    
    if (o->maxIter > 0)
      { /* Perform nonlinear optimization of individual image maps: */
        imsn_optimize_all_pmaps(&(o->image), &mt, &M, o->perspective, o->maxIter, o->maxErr, o->verbose);
        /* Must redo the distortion and scale equalization: */
        if (o->eqDistortion)
          { imsn_equalize_pmap_distortion(&(o->image), &M, o->perspective, o->verbose); }
        imsn_adjust_global_scale(&(o->image), &M, o->minScale, o->maxScale, o->verbose);
      }
   
    /* Adjust rotation, if reqested: */
    if (o->eqRotation)
      { imsn_equalize_rotation(&(o->image), &M, o->verbose); }
    
    /* Adjust translation, if reqested: */
    imsn_adjust_translation(&(o->image), &M, &(o->center), &(o->minCoords), o->verbose); 
                                        
    /* Output and checking: */
    if (o->verbose)
      { imsn_show_all_pmaps(&(o->image), &M, "final");
        imsn_show_all_matches(&(o->image), &mt, &M, "final");
      }

    imsn_output_all_pmaps(o->outPrefix, &(o->image), &M, o->verbose); 
    imsn_output_all_matches(o->outPrefix, &(o->image), &mt, &M, o->verbose); 
    return 0;
  }

void imsn_read_all_data
  ( options_image_vec_t *image, /* (IN) Image options. */
    options_match_vec_t *match, /* (IN) Point pair options. */
    bool_t verbose,             /* (IN) Print diagnostics. */
    match_data_vec_t *mt        /* (OUT) Point pair lists. */  
  )
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    int nmt = match->ne; /* Number of point pair lists. */
    
    /* (Re)allocate the output vector: */
    match_data_vec_trim(mt, nmt);
    
    /* Read the point pair lists: */
    int k;
    for (k = 0; k < nmt; k++)
      { options_match_t *matchk = &(match->e[k]);
        int ixA = imsn_lookup_image_id(matchk->idA, image);
        int ixB = imsn_lookup_image_id(matchk->idB, image);
        r2_vec_t pA = r2_vec_new(0);    /* Reference points in image A. */
        r2_vec_t pB = r2_vec_new(0);    /* Reference points in image B. */
        imsn_read_point_pair_list(matchk->fname, &pA, &pB, verbose);
        assert(pA.ne == pB.ne);
        bool_t AOK = imsn_check_point_containment(&(image->e[ixA]), &pA);
        bool_t BOK = imsn_check_point_containment(&(image->e[ixB]), &pB);
        demand((AOK & BOK), "error(s) in matched point list -- aborted");
        mt->e[k] = (match_data_t) { .ixA = ixA, .pA = pA, .ixB = ixB, .pB = pB };
      }
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }
  
void imsn_compute_initial_affine_pmaps
  ( options_image_vec_t *image,   /* (IN) Image options. */      
    match_data_vec_t *mt,         /* (IN) Point pair options. */ 
    bool_t verbose,               /* (IN) Print diagnostics. */  
    hr2_pmap_vec_t *M             /* (OUT) Affine maps. */  
  )
  {  
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    int nim = image->ne; /* Number of images. */

    /* (Re)allocate the output vector: */
    hr2_pmap_vec_trim(M, nim);
    
    /* The map for image 0 is the identity: */
    r3x3_ident(&(M->e[0].dir));
    r3x3_ident(&(M->e[0].inv));
    
    int i;
    for (i = 1; i < nim; i++)
      { /* Choose for image {i} the map that best matches the previous ones: */
        r2_vec_t pA = r2_vec_new(0);    /* Reference points in image A. */
        r2_vec_t pB = r2_vec_new(0);    /* Reference points in image B. */
        imsn_collect_previous_matches(i, M, mt, &pA, &pB, verbose);
        assert(pA.ne == pB.ne);
        demand(pA.ne >= 3, "poorly connected or badly sorted images");
        M->e[i] = imsn_compute_affine_pmap(&pA, &pB, verbose);
      }
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }

void imsn_equalize_pmap_distortion
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M,
    bool_t perspective,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    fprintf(stderr, "!!! %s not implemented !!!\n", __FUNCTION__);
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }

void imsn_adjust_global_scale
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M,
    double minScale,
    double maxScale,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    // bool_t adjust_scale = 

    //   ( (o->minScale > 0) && (isfinite(o->minScale)) ) ||
    //   ( (o->maxScale > 0) && (isfinite(o->maxScale)) );
    
    fprintf(stderr, "!!! %s not implemented !!!\n", __FUNCTION__);
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }

void imsn_optimize_all_pmaps
  ( options_image_vec_t *image,
    match_data_vec_t *mt,
    hr2_pmap_vec_t *M, 
    bool_t perspective, 
    int maxIter,
    double maxErr, 
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    //  bool_t debug = TRUE;
    //  
    //  int np = pA->ne;
    //  assert(np == pB->ne);
    //  
    //  if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    //  if (np < 4)
    //    { fprintf(stderr, "too few point pairs (%d) for projective map, using affine map\n", np); 
    //      return (*M0);
    //    }
    //  
    //  int nx = 8; /* Number of packed parameters. */
    //  
    //  auto hr2_pmap_t compute_S_map(void);
    //    /* Computes the {S} map used by {unpack_parameters}. */
    //       
    //  hr2_pmap_t S = compute_S_map();  /* Pmap from imageA's domain to canonical square. */
    //  
    //  auto void pack_parameters(hr2_pmap_t *M, int nx, double x[]);
    //    /* Packs the homogeneous projective map {M} to the parameter
    //      vector {x[0..nx-1]}, as the images of the four corners 
    //      of image A. Requires {nx==8}. */
    //       
    //  auto void unpack_parameters(int nx, double x[], hr2_pmap_t *M);
    //    /* Unpacks the parameter vector {x[0..nx-1]} to the 
    //      homogeneous map matrix {M}.  The matrix
    //      will have {M[0,0] == 1}. Requires {nx==8}. */
    //       
    //  auto double goalf(int nx, double x[]);
    //    /* Computes the mean squared distance between the positions of the
    //       mapped points {pA*M.dir} and {pB*M.inv}, given the packed
    //       parameters {x[0..nx-1]}. */
    //  
    //  auto bool_t is_ok(int nx, double x[], double Fx);
    //    /* Returns true if the squared error {Fx} is small enough. */
    //  
    //  double esq = imsn_mean_err_sqr(pA, &(M0->dir), pB, &(M0->inv));
    //  double dMax = 20*sqrt(esq);      /* Search radius around initial guess. */
    //  double dBox = FALSE;             /* ??? Domain is ball, not box. */
    //  double rIni = 0.250*dMax;        /* Initial probe radius. */
    //  double rMin = 0.5;               /* Minimum probe radius. */
    //  double rMax = 0.500*dMax;        /* Maximum probe radius. */
    //  sign_t dir = -1;                 /* Look for minimum. */
    //  
    //  if (verbose) 
    //    { fprintf(stderr, "initial mean error squared = %22.16e\n", esq);
    //      fprintf(stderr, "estimated distance from optimum = %13.6f\n", dMax);
    //      fprintf(stderr, "probe radius = %13.6f [ %13.6f _ %13.6f ]\n", rIni, rMin, rMax);
    //    }
    //    
    //  double x[nx]; /* Initial guess and final optimum parameters. */ 
    //  pack_parameters(M0, nx, x);
    //  if (verbose) 
    //    { fprintf(stderr, "plotting the goal function\n");
    //      FILE *wr = open_write("out/plot.txt", TRUE);
    //      minn_plot_1D_gnuplot(wr, nx, goalf, x, 20, rIni);
    //      fclose(wr);
    //    }
    //  if (verbose) { fprintf(stderr, "optimizing\n"); }
    //  double Fx = goalf(nx, x);
    //  if (verbose) { fprintf(stderr, "initial rms error = %13.6f\n", Fx); }
    //  sve_minn_iterate(nx, &goalf, &is_ok, x, &Fx, dir, dMax, dBox, rIni, rMin, rMax, stop, maxIter, verbose);
    //  if (verbose) { fprintf(stderr, "final rms error = %13.6f\n", Fx); }
    //  
    //  /* Unpack projective map: */
    //  hr2_pmap_t M; /* Optimized map. */
    //  unpack_parameters(nx, x, &M);
    //  if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
    //  return M;
    //  
    //  /* --- internal procs ----------------------------------------------------- */
    //  
    //  hr2_pmap_t compute_S_map(void)
    //    { hr2_point_t hp[4];
    //      int ix, iy;
    //      for (ix = 0; ix < 2; ix++)
    //        { double px = (ix == 0 ? LA->c[0] : HA->c[0]);
    //          for (iy = 0; iy < 2; iy++)
    //            { double py = (iy == 0 ? LA->c[1] : HA->c[1]);
    //              hp[2*iy+ix] = (hr2_point_t){{{ 1, px, py }}};
    //            }
    //        }
    //      hr2_pmap_t S = hr2_pmap_from_four_points(&(hp[0]), &(hp[1]), &(hp[2]), &(hp[3]));
    //      return hr2_pmap_inv(&S);
    //    }
    //  
    //  void pack_parameters(hr2_pmap_t *M, int nx, double x[])
    //    { assert(nx == 8);
    //      int ix, iy;
    //      int k = 0;
    //      for (ix = 0; ix < 2; ix++)
    //        { double px = (ix == 0 ? LA->c[0] : HA->c[0]);
    //          for (iy = 0; iy < 2; iy++)
    //            { double py = (iy == 0 ? LA->c[1] : HA->c[1]);
    //              r2_t p = (r2_t){{ px, py }};
    //              r2_t q; imsn_map_point(&p, &(M->dir), &q);
    //              x[k] = q.c[0]; k++;
    //              x[k] = q.c[1]; k++;
    //            }
    //        }
    //    }
    //  
    //  void unpack_parameters(int nx, double x[], hr2_pmap_t *M)
    //    { assert(nx == 8);
    //      int ix, iy;
    //      int k = 0;
    //      hr2_point_t hp[4];
    //      for (ix = 0; ix < 2; ix++)
    //        { for (iy = 0; iy < 2; iy++)
    //            { double qx = x[k]; k++;
    //              double qy = x[k]; k++;
    //              hp[2*iy + ix] = (hr2_point_t){{{ 1, qx, qy }}};
    //            }
    //        }
    //      hr2_pmap_t R = hr2_pmap_from_four_points(&(hp[0]), &(hp[1]), &(hp[2]), &(hp[3]));
    //      (*M) = hr2_pmap_compose(&S, &R); 
    //    }
    //  
    //  double goalf(int nx, double x[])
    //    {
    //      if (debug) { rn_gen_print(stderr, nx, x, "%8.1f", "\n    [ ", " ", " ]\n"); }
    //      hr2_pmap_t M;
    //      unpack_parameters(nx, x, &M);
    //      /* Mean square error on point lists: */
    //      double esq = imsn_mean_err_sqr(pA, &(M.dir), pB, &(M.inv));
    //      /* Deformation penalty: */
    //      double dsqA = imsn_deform_sqr(LA, HA, &(M.dir));
    //      double dsq2 = imsn_deform_sqr(LB, HB, &(M.inv));
    //      double F = esq + 0.0001*(dsqA + dsqB);
    //      if (debug) 
    //        { fprintf(stderr, "    mean squared error = %13.6e", esq); 
    //          fprintf(stderr, " squared deform = %13.6e %13.6e", dsqA, dsqB);
    //          fprintf(stderr, " function = %13.6e\n", F);
    //        }
    //      return F;
    //      
    //    }
    //    
    //  bool_t is_ok(int nx, double x[], double Fx)
    //    {
    //      return Fx < maxErr*maxErr;
    //    }

    fprintf(stderr, "!!! %s not implemented !!!\n", __FUNCTION__);
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }

void imsn_equalize_rotation
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M, 
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    fprintf(stderr, "!!! %s not implemented !!!\n", __FUNCTION__);
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }

void imsn_adjust_translation
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M,
    r2_t *center,
    r2_t *minCoords,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    bool_t set_ctr = r2_is_finite(center);
    bool_t set_min = r2_is_finite(minCoords);
    demand(!(set_ctr & set_min), "centering and min-coords are mutually exclusive");

    if (set_ctr || set_min)
      { 
        /* Get the bounding box of mapped domains: */
        interval_t B[2];
        imsn_all_mapped_domains_bbox(image, M, B);
        /* Compute the translatio  vector {vec}: */
        r2_t vec; 
        if (set_ctr) 
          { r2_t ctr = (r2_t){{ interval_mid(&(B[0])), interval_mid(&(B[1])) }};
            r2_sub(center, &ctr, &vec);
          }
        else if (set_min) 
          { r2_t lo = (r2_t){{ B[0].end[0], B[1].end[0] }};
            r2_sub(minCoords, &lo, &vec);
          }
        else
          { assert(FALSE); }
        /* Apply the translation: */
        imsn_translate_all_images(M, &vec);
      }
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }

void imsn_all_mapped_domains_bbox(options_image_vec_t *image, hr2_pmap_vec_t *M, interval_t B[])
  { 
    int nim = M->ne;  /* Number of images. */
    
    /* Initialize the intervals with {[+oo _ -oo]}: */
    B[0] = (interval_t){{ +INF, -INF }};
    B[1] = (interval_t){{ +INF, -INF }};
    
    /* Enclose all mapped domain corners */
    int i;
    for (i = 0; i < nim; i++)
      { options_image_t *imi = &(image->e[i]);
        hr2_pmap_t *Mi = &(M->e[i]);
        interval_t Bi[2];
        imsn_mapped_domain_bbox(imi, Mi, Bi);
        B[0] = interval_join(&(B[0]), &(Bi[0]));
        B[1] = interval_join(&(B[1]), &(Bi[1]));
      }
  }

void imsn_mapped_domain_bbox(options_image_t *im, hr2_pmap_t *M, interval_t B[])
  {
    B[0] = (interval_t){{ +INF, -INF }};
    B[1] = (interval_t){{ +INF, -INF }};
    int j;
    for (j = 0; j < 4; j++)
      { r2_t *Pj = &(im->P[j]);
        r2_t Qj; imsn_map_point(Pj, &(M->dir), &Qj);

        B[0].end[0] = fmin(B[0].end[0], Qj.c[0]);
        B[0].end[1] = fmax(B[0].end[1], Qj.c[0]);

        B[1].end[0] = fmin(B[1].end[0], Qj.c[1]);
        B[1].end[1] = fmax(B[1].end[1], Qj.c[1]);
      }
  }

void imsn_translate_all_images(hr2_pmap_vec_t *M, r2_t *vec)
  {
    int nim = M->ne;  /* Number of images. */
    
    hr2_pmap_t T = hr2_pmap_translation(vec);
    int i;
    for (i = 0; i < nim; i++)
      { hr2_pmap_t *Mi = &(M->e[i]);
        (*Mi) = hr2_pmap_compose(Mi, &T);
      }
  }

void imsn_show_all_pmaps
  ( options_image_vec_t *image,
    hr2_pmap_vec_t *M, 
    char *descr
  )
  {
    int nim = image->ne; /* Number of images. */
    
    int i;
    for (i = 0; i < nim; i++)
      { options_image_t *imagei = &(image->e[i]);
        hr2_pmap_t *Mi = &(M->e[i]);
        imsn_show_pmap(Mi, imagei->id, descr);
        imsn_show_mapped_domain(imagei, Mi);
        
        /* !!! Show the min, max scale at corners, and avg scale at center. !!! */

        interval_t Bi[2]; 
        imsn_mapped_domain_bbox(imagei, Mi, Bi);
        imsn_show_bbox_and_center(Bi, "mapped domain");
      }
    /* Print the boundinng box of all images: */ 
    interval_t B[2]; 
    imsn_all_mapped_domains_bbox(image, M, B);
    imsn_show_bbox_and_center(B, "all mapped domains");
  }

void imsn_show_all_matches
  ( options_image_vec_t *image,
    match_data_vec_t *mt, 
    hr2_pmap_vec_t *M, 
    char *descr
  )
  {
    int nim = image->ne; /* Number of images. */
    int nmt = mt->ne;    /* Number of point pair lists. */
    
    int k;
    for (k = 0; k < nmt; k++)
      { match_data_t *mtk = &(mt->e[k]);
        
        int ixA = mtk->ixA;
        r2_vec_t *pA = &(mtk->pA);
        assert((ixA >= 0) && (ixA < nim)); 
        char *idA = image->e[ixA].id;
        r3x3_t *MA = &(M->e[ixA].dir);
        
        int ixB = mtk->ixB;
        r2_vec_t *pB = &(mtk->pB);
        assert((ixB >= 0) && (ixB < nim)); 
        char *idB = image->e[ixB].id;
        r3x3_t *MB = &(M->e[ixB].dir);
        
        imsn_show_match(pA, MA, idA, pB, MB, idB, descr);
      }
  }

void imsn_output_all_pmaps
  ( char *outPrefix,
    options_image_vec_t *image,
    hr2_pmap_vec_t *M, 
    bool_t verbose
  )
  {
    int nim = image->ne; /* Number of images. */
    int i;
    for (i = 0; i < nim; i++)
      { options_image_t *imagei = &(image->e[i]);
        hr2_pmap_t *Mi = &(M->e[i]);
        imsn_output_pmap(outPrefix, imagei->id, Mi, verbose);
      }
  }

void imsn_output_all_matches
  ( char *outPrefix,
    options_image_vec_t *image,
    match_data_vec_t *mt, 
    hr2_pmap_vec_t *M, 
    bool_t verbose
  )
  {
    int nim = image->ne; /* Number of images. */
    int nmt = mt->ne;    /* Number of point pair lists. */
    
    int k;
    for (k = 0; k < nmt; k++)
      { match_data_t *mtk = &(mt->e[k]);
        
        int ixA = mtk->ixA;
        r2_vec_t *pA = &(mtk->pA);
        assert((ixA >= 0) && (ixA < nim)); 
        char *idA = image->e[ixA].id;
        r3x3_t *MA = &(M->e[ixA].dir);
        
        int ixB = mtk->ixB;
        r2_vec_t *pB = &(mtk->pB);
        assert((ixB >= 0) && (ixB < nim)); 
        char *idB = image->e[ixB].id;
        r3x3_t *MB = &(M->e[ixB].dir);
        
        imsn_output_match(outPrefix, pA, MA, idA, pB, MB, idB, verbose);
      }
  }

/* ---------------------------------------------------------------------- */

int imsn_lookup_image_id(char *id, options_image_vec_t *image)
  {
    int nim = image->ne; /* Number of images. */
    int i;
    for (i = 0; i < nim; i++)
      { options_image_t *imagei = &(image->e[i]);
        if (strcmp(id, imagei->id) == 0) { return i; }
      }
    fprintf(stderr, "** image id \"%s\" not defined\n", id);
    demand(FALSE, "aborted");
  }

bool_t imsn_check_point_containment(options_image_t *im, r2_vec_t *p)
  {
    int k;
    /* Compute sides {L[0..3]} of domain: */
    hr2_line_t L[4]; /* Domain sides. */
    for (k = 0; k < 4; k++)
      { hr2_point_t a = hr2_from_r2(&(im->P[k]));
        hr2_point_t b = hr2_from_r2(&(im->P[(k + 1) % 4]));
        L[k] = hr2_join(&a, &b);
      }
    /* Check points for containment: */
    int nbugs = 0;
    int j;
    for (j = 0; j < p->ne; j++)
      { r2_t *pj = &(p->e[j]);
        hr2_point_t hj = hr2_from_r2(pj);
        bool_t inside = TRUE; /* For now. */
        for (k = 0; k < 4; k++)
          { hr2_line_t *Lk = &(L[k]);
            if (hr2_side(&hj, Lk) <= 0) { inside = FALSE; }
          }
        if (! inside)
          { fprintf(stderr, "point %d = ", j);
            r2_gen_print(stderr, pj, "%8.1f", "( ", " ", " )");
            fprintf(stderr, " is not inside the domain of image %s\n", im->id);
            nbugs++;
          }
      }
    return (nbugs == 0);
  }


void imsn_collect_previous_matches
  ( int m, 
    hr2_pmap_vec_t *M, 
    match_data_vec_t *mt, 
    r2_vec_t *pA, 
    r2_vec_t *pB,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    int nim = M->ne;  /* Number of images. */
    int nmt = mt->ne; /* Number of point pair lists. */
    
    if (verbose) { fprintf(stderr, "prev images = %d  matched point lists = %d\n", m, nmt); }
    int np = 0;
    
    int k;
    for (k = 0; k < nmt; k++)
      { match_data_t *mtk = &(mt->e[k]);
        
        int ixA = mtk->ixA;
        assert((ixA >= 0) && (ixA < nim)); 
        r2_vec_t *pkA = &(mtk->pA);
        
        int ixB = mtk->ixB;
        assert((ixB >= 0) && (ixB < nim)); 
        r2_vec_t *pkB = &(mtk->pB);
        
        int ixV = -1;
        r2_vec_t *pkV = NULL;
        r2_vec_t *pkM = NULL;
        if ((ixA < m) && (ixB == m))
          { 
            ixV = ixA; pkV = pkA; pkM = pkB;
          }
        else if ((ixB < m) && (ixA == m))
          { 
            ixV = ixB; pkV = pkB; pkM = pkA;
          }

        if (ixV >= 0)
          { assert(pkV->ne == pkM->ne);
            int npk = pkV->ne; /* Number of point pairs in this list. */
            if (verbose) { fprintf(stderr, "gathering %d pairs from %d to %d\n", npk, m, ixV); }
            hr2_pmap_t *Mk = &(M->e[ixV]);
            r2_vec_expand(pA, np + npk - 1);
            r2_vec_expand(pB, np + npk - 1);
            int i;
            for(i = 0; i < npk; i++)
              { imsn_map_point(&(pkV->e[i]), &(Mk->dir), &(pA->e[np]));
                pB->e[np] = pkM->e[i];
                np++;
              }
          }
      }
    r2_vec_trim(pA, np);
    r2_vec_trim(pB, np);
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }

void imsn_read_point_pair_list(char *fname, r2_vec_t *pA, r2_vec_t *pB, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }
    FILE *rd = open_read(fname, verbose);
    (*pA) = r2_vec_new(20);
    (*pB) = r2_vec_new(20);
    int np = 0; /* Number of data pairs. */

    while (TRUE)
      { /* Skip comments, if any: */
        fget_skip_formatting_chars(rd);
        if (fget_skip_and_test_char(rd, '#')) 
          { int c; 
            do { c = fgetc(rd); } while ((c != EOF) && (c != '\n'));
            continue;
          }
        if (fget_skip_and_test_char(rd, EOF)) { break; }
        /* Try to read another point pair: */
        double xA, yA, xB, yB;
        int nscan = fscanf(rd, " ( %lf %lf ) = ( %lf %lf )", &xA, &yA, &xB, &yB);
        if (nscan != 4) 
          { fprintf(stderr, "line %d: bad input format (nscan = %d)\n", np+1, nscan); exit(1); }
        fprintf(stderr, " ( %6.1f %6.1f ) = ( %6.1f %6.1f )\n", xA,yA,xB,yB);
        r2_vec_expand(pA, np); pA->e[np] = (r2_t){{xA, yA}}; 
        r2_vec_expand(pB, np); pB->e[np] = (r2_t){{xB, yB}}; 
        np++;
      }
    fclose(rd);
    r2_vec_trim(pA, np); 
    r2_vec_trim(pB, np);
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
  }

hr2_pmap_t imsn_compute_affine_pmap(r2_vec_t *pA, r2_vec_t *pB, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- enter %s ---\n", __FUNCTION__); }

    assert(pA->ne == pB->ne);
    int np = pA->ne;
    if (np < 3) { fprintf(stderr, "too few point pairs (%d) for affine mapping\n", np); exit(1); }

    r2_aff_map_t A = r2_aff_map_from_point_pairs(pB, pA);
    hr2_pmap_t M = hr2_pmap_from_mat_and_disp(&(A.mat), &(A.disp));
    if (verbose) { fprintf(stderr, "--- exit %s ---\n", __FUNCTION__); }
    return M;
  }
        
void imsn_map_point(r2_t *p, r3x3_t *M, r2_t *q)
  {
    r3_t hp = (r3_t){{ 1, p->c[0], p->c[1] }};
    r3_t hq;
    r3x3_map_row(&hp, M, &hq);
    double w = hq.c[0];
    double m = fmax(fabs(hq.c[1]), fabs(hq.c[2]));
    if (w <= m*1e-200) 
      { (*q) =  (r2_t){{ INF, INF }}; }
    else
      { (*q) =  (r2_t){{ hq.c[1]/w, hq.c[2]/w }}; } 
  }
    

double imsn_mean_err_sqr(r2_vec_t *pA, r3x3_t *MA, r2_vec_t *pB, r3x3_t *MB)
  {
    int np = pA->ne;
    assert(np == pB->ne);
    
    int k;
    double sum2 = 0.0;
    for (k = 0; k < np; k++)
      {
        r2_t *pAk = &(pA->e[k]);
        r2_t qAk; imsn_map_point(pAk, MA, &qAk);
        r2_t *pBk = &(pB->e[k]);
        r2_t qBk; imsn_map_point(pBk, MB, &qBk);
        double d2 = r2_dist_sqr(&qAk, &qBk);
        sum2 += d2;
      }
    return sum2/np;
  }

double imsn_deform_sqr(r2_t P[], r3x3_t *M)
  {
    /* Compute the logs {logr[0..5]} of ratios of distances between image corners: */
    /* Compute the logs {logr[0..5]} of ratios of distances between image corners: */
    double logr[6];
    int i, j;
    int k = 0;
    for (i = 0; i < 4; i++)
      { for (j = 0; j < i; j++)
          { r2_t uo = P[i];
            r2_t vo = P[j];
            /* Get the endpoints {um,vm} of mapped side/diagonal {k}: */
            r2_t um; imsn_map_point(&uo, M, &um);
            r2_t vm; imsn_map_point(&vo, M, &vm);
            /* Compute and save the log of the distance ratio: */
            double dosq = r2_dist_sqr(&uo, &vo);
            double dmsq = r2_dist_sqr(&um, &vm);
            logr[k] = log(dmsq/dosq);
            k++;
          }
      }
    assert(k == 6);
      
    /* Compute the variance of the logs: */
    double sum = 0;
    for (k = 0; k < 6; k++) { sum += logr[k]; }
    double avg = sum/6;
    double sum2 = 0;
    for (k = 0; k < 6; k++) { double dk = logr[k] - avg; sum2 += dk*dk; }
    double var = sum2/5;
    return var;
  }

void imsn_show_pmap(hr2_pmap_t *M, char *id, char *descr)
  {
    fprintf(stderr, "  image %s - %s matrix (direct):\n", id, descr);
    r3x3_gen_print(stderr, &(M->dir), "%13.6e", "", "\n", "\n", "    [ ", " ", " ]");

    fprintf(stderr, "  image %s - %s matrix (inverse):\n", id, descr);
    r3x3_gen_print(stderr, &(M->inv), "%13.6e", "", "\n", "\n", "    [ ", " ", " ]");
  }

void imsn_show_mapped_domain(options_image_t *im, hr2_pmap_t *M)
  {
    r2_t Q[4]; /* The mapped domain corners: */
    fprintf(stderr, "mapped domain corners = ");
    int j;
    for (j = 0; j < 4; j++)
      { r2_t *Pj = &(im->P[j]);
        r2_t *Qj = &(Q[j]);
        imsn_map_point(Pj, &(M->dir), Qj);
        r2_gen_print(stderr, Qj, "%8.1f", "( ", " ", " )");
      }
    r2_t cod = imsn_domain_centroid(Q);
    fprintf(stderr, " centroid = ");
    r2_gen_print(stderr, &cod, "%8.1f", "( ", " ", " )");       
    fprintf(stderr, "\n");
  }
  
r2_t imsn_domain_centroid(r2_t Q[])
  {
    hr2_line_t L[2];
    int j;
    for (j = 0; j < 2; j++)
      { hr2_point_t a = hr2_from_r2(&Q[j]);
        hr2_point_t b = hr2_from_r2(&Q[j+2]);
        L[j] = hr2_join(&a, &b);
      }
    hr2_point_t u = hr2_meet(&(L[0]), &(L[1]));
    assert(u.c.c[0] != 0);
    r2_t cod = r2_from_hr2(&u);
    return cod;
  }

void imsn_show_bbox_and_center(interval_t B[], char *what)
  {
    fprintf(stderr, "bounding box of %s = ", what);
    interval_gen_print(stderr, &(B[0]), "%8.1f", "[ ", " _ ", " ]");
    fprintf(stderr, " × ");
    interval_gen_print(stderr, &(B[1]), "%8.1f", "[ ", " _ ", " ]");
    r2_t ctr = (r2_t){{ interval_mid(&(B[0])), interval_mid(&(B[1])) }};
    fprintf(stderr, " center = ");
    r2_gen_print(stderr, &ctr, "%8.1f", "( ", " ", " )");       
    fprintf(stderr, "\n");
  }
    
void imsn_output_pmap(char *outPrefix, char *id, hr2_pmap_t *M, bool_t verbose)
  {
    imsn_output_matrix(outPrefix, id, "dir", &(M->dir), verbose);
    imsn_output_matrix(outPrefix, id, "inv", &(M->inv), verbose);
  }
    
void imsn_output_matrix(char *outPrefix, char *id, char *dir, r3x3_t *M, bool_t verbose)
  {
    char *fname = NULL;
    asprintf(&fname, "%s-matrix-%s-%s.txt", outPrefix, id, dir);
    FILE *wr = open_write(fname, verbose);
    fprintf(wr, "# Last edited on DATE TIME by USER\n");
    fprintf(wr, "# Created by %s %s\n", PROG_NAME, PROG_VERS);
    fprintf(wr, "\n");
    r3x3_gen_print(wr, M, "%24.16e", "", "", "", "", " ", "\n");
    fclose(wr);
    free(fname);
  }

void imsn_show_match
  ( r2_vec_t *pA, 
    r3x3_t *MA, 
    char *idA, 
    r2_vec_t *pB, 
    r3x3_t *MB, 
    char *idB, 
    char *descr
  )
  {
    int np = pA->ne;
    assert(np == pB->ne);
    
    fprintf(stderr, "matched points between images %s and %s\n", idA, idB);
    int k;
    double sum2 = 0;
    double emax = -INF;
    for(k = 0; k < np; k++)
      { r2_t *pAk = &(pA->e[k]);
        r2_t qAk; imsn_map_point(pAk, MA, &qAk);
        r2_t *pBk = &(pB->e[k]);
        r2_t qBk; imsn_map_point(pBk, MB, &qBk);
        r2_gen_print(stderr, pAk, "%8.2f", "( ", " ", " )");
        fprintf(stderr, " -> ");
        r2_gen_print(stderr, &qAk, "%8.2f", "( ", " ", " )");
        fprintf(stderr, "  ");
        r2_t d; r2_sub(&qAk, &qBk, &d);
        r2_gen_print(stderr, &d, "%+6.2f", "( ", " ", " )");
        double e2 = r2_dist_sqr(&qAk, &qBk);
        double e = sqrt(e2);
        fprintf(stderr, "%6.2f", e);
        fprintf(stderr, "  ");
        r2_gen_print(stderr, &qBk, "%8.2f", "( ", " ", " )");
        fprintf(stderr, " <- ");
        r2_gen_print(stderr, pBk, "%8.2f", "( ", " ", " )");
        fprintf(stderr, "\n");
        
        sum2 += e2;
        emax = fmax(emax, e);
      }
    double erms = sqrt(sum2/np);
    fprintf(stderr, "root mean squared error = %8.4f\n", erms);
    fprintf(stderr, "max error               = %8.4f\n", emax);
  }

void imsn_output_match
  ( char *outPrefix,
    r2_vec_t *pA, 
    r3x3_t *MA, 
    char *idA, 
    r2_vec_t *pB, 
    r3x3_t *MB, 
    char *idB,
    bool_t verbose
  )
  {
    int np = pA->ne;
    assert(np == pB->ne);
    
    char *fname = NULL;
    asprintf(&fname, "%s-match-%s-%s.txt", outPrefix, idA, idB);
    FILE *wr = open_write(fname, verbose);
    fprintf(wr, "# Last edited on DATE TIME by USER\n");
    fprintf(wr, "# Created by %s %s\n", PROG_NAME, PROG_VERS);
    fprintf(wr, "\n");
    int k;
    for(k = 0; k < np; k++)
      { r2_t *pAk = &(pA->e[k]);
        r2_t qAk; imsn_map_point(pAk, MA, &qAk);
        r2_t *pBk = &(pB->e[k]);
        r2_t qBk; imsn_map_point(pBk, MB, &qBk);

        r2_gen_print(wr, &qAk, "%24.16e", " ", " ", " ");
        r2_gen_print(wr, &qBk, "%24.16e", " ", " ", "\n");
      }
    fclose(wr);
    free(fname);
  }
  
options_t *imsn_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    /* Parse keyword-based arguments: */
    
    /* Image data: */
    int nim = 0;
    while (argparser_keyword_present(pp, "-image"))
      { options_image_t oim;
        oim.id = argparser_get_next_non_keyword(pp);
        hr2_point_t H[4]; /* Corners in homogeneous coordinates. */
        int k;
        for (k = 0; k < 4; k++)
          { r2_t *Pk = &(oim.P[k]);
            Pk->c[0] = argparser_get_next_double(pp, -100000, +100000);
            Pk->c[1] = argparser_get_next_double(pp, -100000, +100000);
            H[k] = hr2_from_r2(Pk);
          }
        /* Check domain for convexity: */
        for (k = 0; k < 4; k++)
          { hr2_point_t *p = &(H[k]);
            hr2_point_t *q = &(H[(k+1)%4]);
            int j;
            for (j = 0; j < 2; j++)
              { hr2_point_t *u = &(H[(k+2+j)%4]);
                if (hr2_orient(p, q, u) <= 0) 
                  { argparser_error(pp, "domain is not ccw convex"); }
              }
          }
        options_image_vec_expand(&(o->image), nim);
        o->image.e[nim] = oim;
        nim++;
      }
    options_image_vec_trim(&(o->image), nim);
                                             
    /* Image pair data: */
    int nmt = 0;
    while (argparser_keyword_present(pp, "-match"))
      { options_match_t omt;
        omt.idA = argparser_get_next_non_keyword(pp);
        omt.idB = argparser_get_next_non_keyword(pp);
        omt.fname = argparser_get_next_non_keyword(pp);
        options_match_vec_expand(&(o->match), nmt);
        o->match.e[nmt] = omt;
        nmt++;
      }
    options_match_vec_trim(&(o->match), nmt);

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    if (argparser_keyword_present(pp, "-perspective"))
      { o->perspective = argparser_get_next_bool(pp);  }
    else
      { o->perspective = TRUE; }

    if (argparser_keyword_present(pp, "-eqDistortion"))
      { o->eqDistortion = argparser_get_next_bool(pp);  }
    else
      { o->eqDistortion = FALSE; }

    if (argparser_keyword_present(pp, "-minScale"))
      { o->minScale  = argparser_get_next_double(pp, 0, 100000); }
    else
      { o->minScale = 0; }
    
    if (argparser_keyword_present(pp, "-maxScale"))
      { o->maxScale = argparser_get_next_double(pp, o->minScale, 100000); } 
    else
      { o->maxScale = 0; }

    if (argparser_keyword_present(pp, "-eqRotation"))
      { o->eqRotation = argparser_get_next_bool(pp);  }
    else
      { o->eqRotation = FALSE; }

    if (argparser_keyword_present(pp, "-center"))
      { o->center.c[0] = argparser_get_next_double(pp, -INF, +INF);
        o->center.c[1] = argparser_get_next_double(pp, -INF, +INF);
      } 
    else
      { o->center = (r2_t){{ INF, INF }}; }

    if (argparser_keyword_present(pp, "-minCoords"))
      { o->minCoords.c[0] = argparser_get_next_double(pp, -INF, +INF);
        o->minCoords.c[1] = argparser_get_next_double(pp, -INF, +INF);
      } 
    else
      { o->minCoords = (r2_t){{ INF, INF }}; }

    if (argparser_keyword_present(pp, "-maxIter"))
      { o->maxIter = argparser_get_next_int(pp, 0, 100000); }
    else
      { o->maxIter = 5; }

    if (argparser_keyword_present(pp, "-maxErr"))
      { o->maxErr = argparser_get_next_double(pp, 0, +INF); }
    else
      { o->maxErr = 0.0; }

    if (argparser_keyword_present(pp, "-verbose"))
      { o->verbose = argparser_get_next_bool(pp);  }
    else
      { o->verbose = FALSE; }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

vec_typeimpl(options_image_vec_t,options_image_vec,options_image_t);
vec_typeimpl(options_match_vec_t,options_match_vec,options_match_t);
vec_typeimpl(match_data_vec_t,match_data_vec,match_data_t);
vec_typeimpl(hr2_pmap_vec_t,hr2_pmap_vec,hr2_pmap_t);


