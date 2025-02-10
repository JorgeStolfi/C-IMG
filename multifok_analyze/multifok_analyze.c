#define PROG_NAME "multifok_analyze"
#define PROG_DESC "Multi-focus stereo microscopy"
#define PROG_VERS "1.0"

/* Last edited on 2025-01-30 07:40:25 by stolfi */

#define multifok_analyze_C_COPYRIGHT \
    "© 2017 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -imageFormat {IMAGE_FORMAT} \\\n" \
  "    -framePattern {FRAME_PATTERN} \\\n" \
  "    { -frameID {FRAME_ID[f]} }.. \\\n" \
  "    -reference {REFIMG} \\\n" \
  "    -windowSize {NW} \\\n" \
  "    -noise {NOISE} \\\n" \
  "    -threshold {THRESH} \\\n" \
  "    [ -verbose ] \\\n" \
  "    -outDir {OUT_DIR} \\\n"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "   Reads a stack of {NF} images of an object (/frames/), taken with a camera with" \
  " limited depth of focus, all from the same direction and with same" \
  " lighting, aperture, and speed, but with different object-camera" \
  " distances.  Reads also a single image where all parts of the object" \
  " are in focus. Tries to determine the focus detection operator that would give" \
  " the most consistent results for those images.\n" \
  "\n" \
  "   The input files must be in the specified format, all with" \
  " the same size, and must be aligned and scaled so that the" \
  " focused parts of the object appear in the same scale in every image.\n" \
  "\n" \
  "THEORY\n" \
  "   The focus operator is determined by linear regression.  The dependent" \
  " variable {s2[p]} for each observation {p} is the local quadratic similarity" \
  " between the reference image and each frame image.  The independent" \
  " variables are a set of {NQ} /terms/ {Q[p][0..NQ-1]} -- local quadratic" \
  " \"energy\" functionals of the frame images. The idea is that, the larger" \
  " these quadratic functionals are at a given pixel, the more likely that" \
  " the frame at that pixel is in focus -- that is, it matches the perfectly" \
  " focused image.\n" \
  "\n" \
  "   More precisely, the program sweeps a sliding window {NW} by {NW} pixels" \
  " over the frame images and the reference image. For each window position {d}, the program extracts" \
  " from the reference image the {NS=NW*NW} samples in that window as a" \
  " vector {rr[0..NS-1]}.  If the deviation of those samples is less than {THRESH}, that" \
  " window postion is considered not significant, and ignored.  Otherwise the following steps are applied" \
  " to that window position, yielding {NF} observations.\n" \
  "\n" \
  "  Thus, if there are {ND} significant window positions, the number {NP} of" \
  " observations for the linear regression is {ND*NF}.  The observation that is obtained from significant" \
  " window position number {d \\in 0..ND-1} and frame index {f \\in 0..NF-1} gets" \
  " an observation index {p = d*NF + f \\in 0..NP-1}.\n" \
  "\n" \
  "  For each window position {d}, the program normalizes the reference sample vector {rr[0..NS-1]} to zero mean and" \
  " unit norm.  Then, for each frame index {f} in {0..NF-1}, it extracts the" \
  " samples from the same window position {d} in that frame image as a vector" \
  " {fr[0..NB-1]}, and applies to them the same normalization" \
  " applied to {rr[0..NS-1]}.  The /local discrepancy/ {e2[p]} of the observation" \
  " is defined as {D2/4} where {D2} is the Euclidean distance" \
  " squared between the normalized vectors {rr} and {ff}.  The value" \
  " of the independent variable for that window position is" \
  " the /local congruence/ {c2[p] = 1 - e2[p]}.  All these operations are weighted" \
  " by the window weight mask defined in {multifok_focus_op_prod_weights}.\n" \
  "\n" \
  "  The vector {fr[0..NS-1]} is then converted" \
  " to another vector {gr[0..NS-1]} by an orthonormal linar map described" \
  " in {multifok_focus_op_basis}.  Then the program computes" \
  " a certain number {NQ} of quadratic functions" \
  " from the vector {gr[0..NS-1]}.  These terms are the independent" \
  " variables of the observation, namely {Q[p][0..NQ-1]}.\n" \
  "\n" \
  "  The (hopefully) optimal focus operator is obtained by linear" \
  " regression of these observations.  It results in a coefficient" \
  " vector {C[0..NQ-1]} of the linear combination {m2[p] = SUM{C[k]*Q[p][k]}} of the quadratic" \
  " functionals that best fits the quadratic discrepancies {e2[p]}.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The program writes to \"{OUT_DIR}/coeff.txt\" a table of coefficients" \
  " for the best focus operator found.   Each line has the term index {q} in {0..NQ-1}, the" \
  " coefficient value, and a description of the term {Q[][q]}.\n" \
  "\n" \
  "  The program writes to \"{OUT_DIR}/data.txt\" the observations used in the" \
  " linear regression.\n" \
  "  " multifok_analyze_write_observations_INFO "\n" \
  "\n" \
  "  The program also writes, for each input frame, four image" \
  " files \"{OUT_DIR}/{NNNNN}_e2.png\", \"{OUT_DIR}/{NNNNN}_m2.png\"," \
  " \"{OUT_DIR}/{NNNNN}_s2.png\", and  \"{OUT_DIR}/{NNNNN}_s2.png\", showing respectively" \
  " the pixel-by-pixel observed discrepancies {e2[p]} between reference image and each frame, the discrepancies {m2[p]} estimated by regression, the normalized similarities {s2[p]}, and the estimated misfocusing {un[p]}. The" \
  " field {NNNNN} is the numeric frame ID specified with the" \
  " corresponding \"-frameID\" argument.  The values are" \
  " mapped linearly from {[0_1]} to {0..65535}.\n" \
  "\n" \
  "  The program also writes, for each frame index {f \\in 0..NF-1} and" \
  " each index {q \\in 0..NQ-1} of a quadratic" \
  " functional used in the" \
  " regression, an image file \"{OUT_DIR}/{NNNNN}_q{QQ}.png\" showing" \
  " the value of that functional at each pixel.  The field {NNNNN} is the frame ID, as" \
  " above, and {QQ} is the functional index {q} as" \
  " two zero-padded decimal digits.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -framePattern {FRAME_PATTERN}\n" \
  "    This mandatory argument specifies the pattern for frame file" \
  " names.  It must include exacly one unquoted {printf} formatting" \
  " code for 32-bit signed integer, like '%06d'.  That code will be" \
  " replaced by the frame number specified in the \"-frameID\" argument.\n" \
  "\n" \
  "  -frameID {FRAME_ID[f]} \n" \
  "    This argument specifies the numeric ID of a partially defocused" \
  " frame.  This command should be repeated for" \
  " each frame in the stack.  The number is parsed in base 10 and" \
  " should be non-negative.\n" \
  "\n" \
  "  -imageFormat {IMAGE_FORMAT}\n" \
  "    This mandatory argument specifies the format of input" \
  " image files.  " image_file_format_arg_INFO "\n" \
  "\n" \
  "  -reference {REFIMG} \n" \
  "    This argument specifies the filename of the image with" \
  " the object entirely in focus.\n" \
  "\n" \
  "  -noise {NOISE}\n" \
  "    This mandatory argument specifies a noise level to be assumed in the images.\n" \
  "\n" \
  "  -threshold {THRESH}\n" \
  "    This mandatory argument specifies the window significance threshold.  Windows" \
  " of the reference image whose sample deviation is less than {THRESH} will be skipped.\n" \
  "\n" \
  "  -windowSize {NW}\n" \
  "    This argument specifies the operator window size, which must be odd.  Note" \
  " that the memory used is proportional to {NW^4} and the processing time is proportional to {NW^8}.\n" \
  "\n" \
  "  -outDir {OUT_DIR}\n" \
  "    This mandatory argument specifies the directory were all output files will be written to.\n" \
  "\n" \
  "  -verbose\n" \
  "    If present, this option requests debugging printouts to {stderr}.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  If you can't see that ALSO, you may need new glasses.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created Sep/2017 by J.Stolfi.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2017-01-01 Added discrepancy image output. [J.Stolfi]\n" \
  "  2017-01-03 Added \"-framePattern\" option. [J.Stolfi]\n" \
  "  2018-09-05 Added output images with quadratic functionals and fitted indicator. [J.Stolfi]\n" \
  "  2018-09-05 Changed definition of {e2} to be discrepancy rather than similarity. [J.Stolfi]\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " multifok_analyze_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <vec.h>
#include <jsfile.h>
#include <argparser.h>

#include <float_image.h>

#include <float_image_read_gen.h>
#include <float_image_write_gen.h>
#include <sample_conv.h>
#include <image_file_format.h>

#include <multifok_focus_op.h>

#include <multifok_analyze_collect_observations.h>
#include <multifok_analyze_compute_rel_similarities.h>
#include <multifok_analyze_fit_coeffs.h>
#include <multifok_analyze_read_frames.h>
#include <multifok_analyze_read_reference.h>

typedef struct options_t 
  { image_file_format_t imageFormat; /* Format of frame files. */
    char *framePattern;              /* Filename pattern of the input images of the stack.. */
    int32_vec_t frameID;             /* Numerical indices of the input images of the stack. */
    char *reference;                 /* Filename of the perfectly focused image. */
    int32_t windowSize;              /* Window width and height. */
    double noise;                    /* Noise level to assume in image. */
    double threshold;                /* Min sample deviation to consider window significant. */
    char *outDir;                    /* Directory for output files. */
    bool_t verbose;                  /* TRUE to print debugging info. */
  } options_t;

options_t *multifok_analyze_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */

float_image_t **multifok_analyze_alloc_indicator_images(int32_t NF, int32_t NX, int32_t NY);
  /* ALlocates {NF} monchrome images with {NX} columns and {NY} rows of pixels,
    and fills them with zeros. */

void multifok_analyze_write_indicator_images
  ( char *outDir,
    int32_t NF, 
    int32_t frameID[],
    int32_t ND,
    int32_t ix[], 
    int32_t iy[], 
    int32_t NP, 
    double e2[], 
    double m2[], 
    double s2[], 
    double un[], 
    float_image_t *dbimg[],
    bool_t verbose
  );
  /* Writes debugging images with the values of the indicators
    {e2[p],m2[p],s2[p],un[p]} for each of {NF} frames.
    
    Assumes that that {dbimg} is a vector of {NF} pointers to images, all
    monochromatic with the same size as the input frames, with samples
    in {[0_1]}.
    
    Fills each image {dbimg[f]} with the data from {e2,m2,s2,un},
    then writes  to files "{outDir}/{NNNNN}_{tag}.png" where {NNNNN}
    is the numeric frame ID {frameID[f]} formatted as 5 digits, zero-padded, and
    {tag} is "e2", "m2", "s2", and "un", respectively. */

void multifok_analyze_set_debug_images
  ( int32_t NF, 
    int32_t frameID[],
    float_image_t *dbimg[], 
    int32_t ND, 
    int32_t ix[], 
    int32_t iy[], 
    int32_t NP, 
    double ind[],
    bool_t verbose
  );
  /* Assumes that {NF} is the number of input frames, with frame IDs
    {frameID[0..NF-1]}, and that {dbimg} is a vector of {NF} pointers to
    images, all monochromatic with the same size as the input frames.
    
    Assumes that {ND,ix,iy,NP} are as described for 
    {multifok_analyze_collect_observations}, and {ind} is a vector
    of {NP = ND*NF} indicator values such as {e2}, {m2}, {s2}, or {un}.  
    Specifically, for each {d} in {0..ND-1} and each {f} in {NF-1}, stores 
    {ind[p]} in the pixel {ix[d],iy[d]} of the proper {dbimg[f]}; where {p=d*NF + f}. 
    
    If {verbose} is true, prints the min and max values of {ind[p]} for each frame. */

void multifok_analyze_write_debug_images
  ( char *outDir,
    char *tag,
    int32_t NF, 
    int32_t frameID[], 
    float_image_t *dbimg[],
    double v0,
    double vM
  );
  /* Assumes that {NF} is the number of input frames, and that {dbimg}
    is a vector of {NF} pointers to images, all monochromatic with the
    same size as the input frames, with samples in {[0_1]}.
    
    Writes each image {dbimg[f]}  to files "{outDir}/dbimg_{NNNNN}.png" where {NNNNN}
    is the numeric frame ID {frameID[f]} formatted as 5 digits, zero-padded. 
    
    The sample values are normalized so that the values {v0} and {vM}
    map to 1.00 and 0.02, respectively. */

void multifok_analyze_write_quadratic_term_images
  ( char *outDir, 
    int32_t NF, 
    int32_t frameID[], 
    int32_t NQ,
    float_image_t *qtimg[]
  );
  /* Assumes that {NF} is the number of input frames, and that {qtimg}
    is a vector of {NF*NQ} pointers to images, all monochromatic with the
    same size as the input frames, with samples in {[0_1]}.  Specifically,
    assumes that image {qtimg[f*NQ + q]} shows the quadratic factor {q}
    extracted from frame {f}.
    
    Writes each image {qtimg[f*NQ + q]}  to files "{outDir}/qtimg_{NNNNN}_q{QQ}.png" where {NNNNN}
    is the numeric frame ID {frameID[f]}, formatted as 5 digits, and {QQ} is 
    the quadratic functional index {q}, formatted as 2 digits, both zero-padded. */

void multifok_analyze_write_observations
  ( char *outDir, 
    int32_t NF,
    int32_t frameID[], 
    int32_t ND, 
    int32_t ix[], 
    int32_t iy[], 
    int32_t NP, 
    double e2[], 
    double m2[], 
    double s2[],
    double un[],
    int32_t NQ, 
    double terms[]
  );
  /* Assumes that {NF,ND,ix,iy,NP,e2,NQ,terms} are as described for 
    {multifok_analyze_collect_observations}.
    
    Writes to "{outDir}/data.txt" the regression observations
    {e2,terms}. Namely, each line will describe one observation {p}
    where {p=d*NF+f}, for each {d} in {0..ND-1} and each {f} in
    {0..NF-1}. See {multifok_analyze_write_observations_INFO}.  */
    
#define multifok_analyze_write_observations_INFO \
  "The output will consist of {ND} blocks of data lines, one block" \
  " for each of {ND} window positions in the reference image, separated" \
  " by blank lines. Each block will have {NF} consecutive lines, one" \
  " for each of the {NF} input frames.  Thus there will be {NP = ND*NF} data" \
  " lines, each corresponding to one observation used in the linear regression.\n" \
  "\n" \
  "   Each line will have the format\n" \
  "\n" \
  "    \"{p} {d} {ix[d]} {iy[d]} {f} {frID[f]} {e2[p]} {m2[p]} {s2[p]} {un[p]} {terms[p][0..NQ-1]}\"\n" \
  "\n" \
  "  where\n" \
  "\n" \
  "    {d} is the window position index in {0..ND-1},\n" \
  "\n" \
  "    {ix[d]} and {iy[d]} are the column and row of the window's center pixel,\n" \
  "\n" \
  "    {f} is the frame index in {0..NF-1},\n" \
  "\n" \
  "    {frID[f]} is the frame's numeric ID,\n" \
  "\n" \
  "    {p} is the index in {0..NP-1} of a significant observation,\n" \
  "\n" \
  "    {e2[p]} is the normalized quadratic discrepancy between {WF[p]} and {WR[d]},\n" \
  "\n" \
  "    {m2[p]} is the approximation of {e2[p]} estimated from {WF[p]} alone by the fitted formula,\n" \
  "\n" \
  "    {s2[p]} is the similarity derived from {e2[p]} or {m2[p]},\n" \
  "\n" \
  "    {un[p]} is the misfocusing derived from {e2[p]} or {m2[p]}, and\n" \
  "\n" \
  "    {terms[p][0..NQ-1]} are the quadratic functionals {Q[p][0..NQ-1]} obtained from {WF[p]}.\n" \
  "\n" \
  "  Each block corresponds to one window in the reference frame, centered" \
  " at column {ix[d]} and row {iy[d]}, whose samples {WR[d]} are sufficiently" \
  " non-uniform.  The values in each data line are derived from those samples" \
  " and from the samples {WF[p]} of the same window in frame image number {f}.\n" \
  "\n" \
  "  The observed discrepancy {e2[p]} is the squared Euclidean distance {D2[p]} between" \
  " the samples of {WF[p]} and {WR[d]}.\n" \
  "\n" \
  "  The quadratic functionals {Q[p][0..NQ-1]} are extracted from the frame samples {WF[p]}.\n" \
  "\n" \
  "  The estimated discrepancy {m2[p]} is a linear combination of the" \
  " quadratic functionals {Q[p][0..NQ-1]} with the coefficients obtained" \
  " by linear regression of {e2[p]} against those functionals.\n" \
  "\n" \
  "  The estimated similarity {s2[p]} is obtained from either the observed" \
  " discrepancy {e2[p]} or from its estimate {m2[p]}.  The discrepancy" \
  " values are complemented, shifted, and rescaled so that the" \
  " approximate maximum within each data block is mapped to 0 and" \
  " the approximate minimum within that block is mapped to 1, with outliers" \
  " clipped to those values.\n" \
  "\n" \
  "  The estimated misfocusing {un[p]} is obtained by estimating" \
  " the (fractional) index {fopt} such that {s2f[fopt]} is the maximum" \
  " of {s2f[0..N-1]}, and subtracting that from the frame index {f}."
  
void multifok_analyze_write_coeffs(char *outDir, int32_t NQ, char *tname[], double coeff[]);
  /* Writes to "{outDir}/coeff.txt" the coefficients {coeff[0..NQ-1]} 
    computed by {multifok_analyze_fit_coeffs}, with the corresponding descriprions {tname[0..NQ-1]}. */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = multifok_analyze_parse_options(argc, argv);
    
    /* Get the reference image: */
    if (o->verbose) { fprintf(stderr, "reading reference image...\n"); }
    float_image_t *ref = multifok_analyze_read_reference(o->reference, o->imageFormat, o->verbose);

    int32_t NC = (int32_t)ref->sz[0];
    int32_t NX = (int32_t)ref->sz[1];
    int32_t NY = (int32_t)ref->sz[2];
    demand(NC == 1, "images must be monochromatic");
    
    /* Get the images: */
    if (o->verbose) { fprintf(stderr, "computing focus operator tables...\n"); }
    int32_t NF = o->frameID.ne;
    int32_t *frameID = o->frameID.e;
    float_image_t **frame = multifok_analyze_read_frames
      ( o->framePattern, NF, frameID, o->imageFormat, NC, NX, NY, o->verbose);
    
    /* Collect and write out the observations: */
    if (o->verbose) { fprintf(stderr, "collecting the focus operator observations...\n"); }
    int32_t NW = o->windowSize;
    int32_t ND = INT_MIN; /* Number of significant window positions found. */
    int32_t NP = INT_MIN; /* Number of observations collected. */
    int32_t NQ = INT_MIN; /* Number of quadratic terms per observation. */
    char **tname = NULL;  /* Descriptions of the terms, {NQ} elements. */
    int32_t *ix = NULL;   /* The column indices of the significant window centers. */
    int32_t *iy = NULL;   /* The row indices of the significant window centers. */
    double *terms = NULL; /* The matrix quadratic terms, {NP} rows of {NQ} columns. */
    double *e2 = NULL;    /* The vector of quadratic discrepancies to be fitted, {NP} elements. */
    float_image_t **qtimg = NULL; /* Images with quadratic terms. */
    multifok_analyze_collect_observations
      (ref, NF, frame, NW, o->noise, o->verbose, &ND, &ix, &iy, &NP, &e2, &NQ, &tname, &terms, &qtimg);
    fprintf(stderr, "input frames NF = %d\n", NF);
    fprintf(stderr, "window positions ND = %d\n", ND);
    fprintf(stderr, "observations NP = %d\n", NP);
    fprintf(stderr, "quadratic terms NQ = %d\n", NQ);
     
    /* Write the quadratic term images: */
    if (o->verbose) { fprintf(stderr, "writing the quadratic term images...\n"); }
    multifok_analyze_write_quadratic_term_images(o->outDir, NF, frameID, NQ, qtimg);

    /* Compute the optimal coefficients and fitted values: */
    if (o->verbose) { fprintf(stderr, "fitting optimal coefficients...\n"); }
    double *coeff = NULL;  /* Coefficients of the fitted formula, {NQ} elements. */
    double *m2 = NULL;     /* The vector of fitted values, {NP} elements. */
    multifok_analyze_fit_coeffs(NP, NQ, terms, e2, &coeff, &m2); 

    /* Write the fitted coefficients: */
    if (o->verbose) { fprintf(stderr, "writing regression coefficients...\n"); }
    multifok_analyze_write_coeffs(o->outDir, NQ, tname, coeff);
    
    /* Compute the similarities and estimated misfocusing: */
    double *s2 = NULL;    /* The relative similarities, with {NP} elements. */
    double *un = NULL;    /* The estimated misfocus, with {NP} elements. */
    bool_t debug_rel_simil = TRUE; /* TRUE to debug the computation of relative similarities. */
    if (o->verbose) { fprintf(stderr, "computing similarity and misfocusing indicators...\n"); }
    multifok_analyze_compute_rel_similarities(NF, ND, e2, &s2, &un, debug_rel_simil);
    
    /* Write the debugging images: */
    float_image_t **dbimg = multifok_analyze_alloc_indicator_images(NF, NX, NY);
    multifok_analyze_write_indicator_images(o->outDir, NF, frameID, ND, ix, iy, NP, e2, m2, s2, un, dbimg, o->verbose);

    /* Write the observations: */
    if (o->verbose) { fprintf(stderr, "writing the observation data...\n"); }
    multifok_analyze_write_observations(o->outDir, NF, frameID, ND, ix, iy, NP, e2, m2, s2, un, NQ, terms);

    if (o->verbose) { fprintf(stderr, "done.\n"); }
    return 0;
  }

float_image_t **multifok_analyze_alloc_indicator_images(int32_t NF, int32_t NX, int32_t NY)
  { float_image_t **img = notnull(malloc(NF*sizeof(float_image_t *)), "no mem");
    for (uint32_t i = 0;  i < NF; i++)
      { img[i] = float_image_new(1, NX, NY); }
    return img;
  }

void multifok_analyze_write_indicator_images
  ( char *outDir,
    int32_t NF, 
    int32_t frameID[],
    int32_t ND,
    int32_t ix[], 
    int32_t iy[], 
    int32_t NP, 
    double e2[], 
    double m2[], 
    double s2[], 
    double un[], 
    float_image_t *dbimg[],
    bool_t verbose
  )
  {
    demand(NP == ND*NF, "inconsistent counts of window positions, observations, and frames");

    if (verbose) { fprintf(stderr, "writing the observed discrepancy {e2} images...\n"); }
    multifok_analyze_set_debug_images(NF, frameID, dbimg, ND, ix, iy, NP, e2, verbose);
    multifok_analyze_write_debug_images(outDir, "e2", NF, frameID, dbimg, 0.0, 1.0);

    if (verbose) { fprintf(stderr, "writing the estimated discrepancy {m2} images...\n"); }
    multifok_analyze_set_debug_images(NF, frameID, dbimg, ND, ix, iy, NP, m2, verbose);
    multifok_analyze_write_debug_images(outDir, "m2", NF, frameID, dbimg, 0.0, 1.0);

    if (verbose) { fprintf(stderr, "writing the normalized similarity {s2} images...\n"); }
    multifok_analyze_set_debug_images(NF, frameID, dbimg, ND, ix, iy, NP, s2, verbose);
    multifok_analyze_write_debug_images(outDir, "s2", NF, frameID, dbimg, 0.0, 1.0);

    if (verbose) { fprintf(stderr, "writing the misfocusing indicator {un} images...\n"); }
    multifok_analyze_set_debug_images(NF, frameID, dbimg, ND, ix, iy, NP, un, verbose);
    multifok_analyze_write_debug_images(outDir, "un", NF, frameID, dbimg, 0.0, 1.0);
  }

void multifok_analyze_set_debug_images
  ( int32_t NF, 
    int32_t frameID[],
    float_image_t **dbimg, 
    int32_t ND, 
    int32_t ix[], 
    int32_t iy[], 
    int32_t NP, 
    double ind[],
    bool_t verbose
  )
  {
    demand(NP == ND*NF, "inconsistent counts of window positions, observations, and frames");
    
    for (uint32_t f = 0;  f < NF; f++)
      { /* Non-significant positions are set to black: */
        float_image_fill_channel(dbimg[f], 0, 0.0f);
        double vmin = +INF;  /* Min signif value for this frame. */
        double vmax = -INF;  /* Max signif value for this frame. */
        for (uint32_t d = 0;  d < ND; d++)
          { int32_t p = d*NF + f; /* Index of observation. */
            double indp = ind[p];
            if (indp < vmin) { vmin = indp; }
            if (indp > vmax) { vmax = indp; }
            float val = (float)indp;
            if (val > 0.99f) { val = 0.99f; }
            if (val < 0.01f) { val = 0.01f; }
            int32_t ixd = ix[d];
            int32_t iyd = iy[d];
            float_image_set_sample(dbimg[f], 0, ixd, iyd, val);
          }
        if (verbose)
          { fprintf(stderr, "frame %05d values in [ %10.7f _ %10.7f ]\n", frameID[f], vmin, vmax); }
      }
   }

void multifok_analyze_write_debug_images
  ( char *outDir, 
    char *tag,
    int32_t NF, 
    int32_t frameID[],
    float_image_t **dbimg,
    double v0,
    double vM
  )
  {
    for (uint32_t f = 0;  f < NF; f++)
      { char *fname = jsprintf("%s/%05d_%s.png", outDir, frameID[f], tag);
        image_file_format_t ffmt = image_file_format_PNG;
        bool_t yUp = TRUE;
        double gammaEnc = 1.0;
        double bias = 0.0;
        bool_t verbose = FALSE;
        float_image_write_gen_named(fname, dbimg[f], ffmt, yUp, v0, vM, gammaEnc, bias, verbose);
        free(fname);
      }
  }
    
void multifok_analyze_write_quadratic_term_images
  ( char *outDir, 
    int32_t NF, 
    int32_t frameID[], 
    int32_t NQ,
    float_image_t **qtimg
  )
  {
    for (uint32_t f = 0;  f < NF; f++)
      { for (uint32_t q = 0;  q < NQ; q++)
          { char *fname = jsprintf("%s/simil_%05d_q%02d.png", outDir, frameID[f], q);
            image_file_format_t ffmt = image_file_format_PNG;
            bool_t yUp = TRUE;
            double v0 = 0.0;
            double vM = 1.0;
            double gammaEnc = 1.0;
            double bias = 0.0;
            bool_t verbose = FALSE;
            float_image_write_gen_named(fname, qtimg[f*NQ + q], ffmt, yUp, v0, vM, gammaEnc, bias, verbose);
            free(fname);
          }
      }
  }
  
void multifok_analyze_write_observations
  ( char *outDir, 
    int32_t NF, 
    int32_t frameID[],
    int32_t ND, 
    int32_t ix[], 
    int32_t iy[], 
    int32_t NP, 
    double e2[], 
    double m2[], 
    double s2[],
    double un[],
    int32_t NQ, 
    double terms[]
  )
  {
    char *fname = jsprintf("%s/data.txt", outDir);
    FILE *wr = open_write(fname, TRUE);

    demand(NP == ND*NF, "inconsistent counts of window positions, observations, and frames");

    /* Write all observations: */
    for (uint32_t d = 0;  d < ND; d++)
      { /* Insert a blank line before the observations of the same window, for {gnuplot}: */
        if (d > 0) { fprintf(wr, "\n"); }
        /* Write observations of window position {d}: */
        for (uint32_t f = 0;  f < NF; f++)
          { int32_t p = d*NF + f; /* Index of observation. */
            /* Write observation: */
            fprintf(wr, "%8d  %5d %5d", d, ix[d], iy[d]);
            fprintf(wr, "  %5d %05d", f, frameID[f]);
            fprintf(wr, "  %8d", p);
            fprintf(wr, " %18.8f", e2[p]);
            fprintf(wr, " %18.8f", m2[p]);
            fprintf(wr, "  %18.8f", s2[p]);
            fprintf(wr, " %6.2f", un[p]);
            double *ts = &(terms[p*NQ]);
            for (uint32_t q = 0;  q < NQ; q++)
              { fprintf(wr, " %10.6f", ts[q]); }
            fprintf(wr, "\n");
          }
      }
    fclose(wr);
    free(fname);
  }

void multifok_analyze_write_coeffs(char *outDir, int32_t NQ, char *tname[], double coeff[])
  {
    char *fname = jsprintf("%s/coeff.txt", outDir);
    FILE *wr = open_write(fname, TRUE);
    
    /* Write coefficients: */
    for (uint32_t q = 0;  q < NQ; q++)
      { fprintf(wr, "%8d %+18.8f %s\n", q, coeff[q], tname[q]); }
    fclose(wr);
    free(fname);
  }

options_t *multifok_analyze_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    argparser_get_keyword(pp, "-framePattern");
    o->framePattern = argparser_get_next_non_keyword(pp);

    /* Parse keyword parameters: */
    o->frameID = int32_vec_new(100);
    int32_t NI = 0;
    while (argparser_keyword_present(pp, "-frameID"))
      { int32_vec_expand(&(o->frameID), NI);
        o->frameID.e[NI] = (int32_t)argparser_get_next_int(pp, 0, INT_MAX); 
        NI++;
      }
    if (NI == 0) { argparser_error(pp, "must specifiy at least one \"-frameID\""); }
    int32_vec_trim(&(o->frameID), NI);

    o->imageFormat = image_file_format_arg_parse(pp, "-imageFormat");

    argparser_get_keyword(pp, "-reference");
    o->reference = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-windowSize");
    o->windowSize = (int32_t)argparser_get_next_int(pp, 3, 7);
        
    argparser_get_keyword(pp, "-noise");
    o->noise = argparser_get_next_double(pp, 0.000, 1.000);

    argparser_get_keyword(pp, "-threshold");
    o->threshold = argparser_get_next_double(pp, 0.000, 1.000);

    o->verbose = argparser_keyword_present(pp, "-verbose");

    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next(pp);

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }


