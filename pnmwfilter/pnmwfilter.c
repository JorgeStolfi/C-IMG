#define PROG_NAME "pnmwfilter"
#define PROG_DESC "PPM/PGM/PBM window filters -- average,deviation,variance,median,percentile,rank,normal,stretch"
#define PROG_VERS "1.0"

#define pnmwfilter_C_COPYRIGHT \
  "Copyright © 2006 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-03-07 19:34:32 by stolfi */

/* !!! Replace the {RMAG} parameter by -scale and -offset options. !!! */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -filter\\\n" \
  "    { median \\\n" \
  "    | percentile {LEVEL} \\\n" \
  "    | average \\\n" \
  "    | deviation \\\n" \
  "    | variance \\\n" \
  "    | normalize {RMAG} \\\n" \
  "    | stretch {LOLEVEL} {HILEVEL} \\\n" \
  "    | rank {LOLEVEL} {HILEVEL} \\\n" \
  "    } \\\n" \
  "  -weights {WFILE} \\\n" \
  "  [ -excludeSelf ] \\\n" \
  "  [ -mask {MFILE} ] \\\n" \
  "  [ -maxval {OMAXVAL} ] \\\n" \
  "  [ -noise {NOISE} ] \\\n" \
  "  [ -badIn {IBAD} [ -keepBad ] ] \\\n" \
  "  [ -badOut {OBAD} ] \\\n" \
  "  [ -isMaskIn {IMK} ] \\\n" \
  "  [ -isMaskOut {OMK} ] \\\n" \
  "  [ -replicate ] \\\n" \
  "  [ INFILE [ OUTFILE ] ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program takes an input PBM, PGM or PPM image {U} and outputs another" \
  " image {V} with the same number of channels, same dimensions, and same sample" \
  " encoding (`raw' or `plain').  For monochromatic images, each pixel {vs} of {V} is" \
  " computed from the samples of {U} in a window centered at the" \
  " corresponding pixel {us} of {U}.  For color (PPM) images," \
  " each channel is processed independently, as a separate monochromatic image.\n" \
  "\n" \
  "FILTER TYPES\n" \
  "  The filter type is selected by the \"-filter\" option.  The" \
  " program currently implements the popular running median filter, the" \
  " general percentile filter (which includes grayscale erosion and" \
  " dilation), average and deviation of window samples, and" \
  " localized versions of image normalization (based on max-min" \
  " or mean-deviation) and image equalization.\n" \
  "\n" \
  "WINDOW IMAGE FILE\n" \
  "  The window size and shape are specified by a user-given" \
  " PGM or PBM image file, the /window image/, with" \
  " odd width and height.  During processing, the center" \
  " pixel of the window image is aligned with the pixel {p} being" \
  " processed.  When computing window statistics such as" \
  " median, mean, or standard deviation, the program interprets each" \
  " window image sample as a numeric `weight' or `area'" \
  " of the input image pixel that is aligned with it.  Thus, for example, if" \
  " the window sample values are 10, 20, 10, 20, 40, and the input image sample" \
  " values aligned with them are 40, 50, 60, 200, 201, then the" \
  " median value is 200, and the 0.20 percentile is 50.\n" \
  "\n" \
  "  Note that some filters, especially \"percentile\" and \"stretch\", will" \
  " give better results if the window weights decay smoothly to 0" \
  " as one moves from the center to the window's borders.  If the window image" \
  " has step-like discontinuities, including at its borders, these will" \
  " usually lead to similar discontinuities in the output image.\n" \
  "\n" \
  "IMAGE MASK\n" \
  "  Besides the mandatory window image, the program also accepts" \
  " an optional PGM or PBM /mask image/ file, with the same size as the input image.  If this" \
  " mask is given, each sample is interpreted as a second weight permanently" \
  " associated to the corresponding input image pixel.  In that case, the effective pixel weights" \
  " used in the filter computations are the products of these mask weights and of the" \
  " window weights currently aligned with them.\n" \
  "\n" \
  "IMAGE ENCODING\n" \
  "  The samples in the input image, weight image, and output image are assumed to have linear" \
  " encoding (gamma 1.0).  The input image samples are converted to float values" \
  " in {[0_1]}, as determined by the input image's {maxval} parameter and" \
  " the \"-badIn\" and \"-isMaskIn\" command line arguments.  Likewise," \
  " the rounding of the computed samples for output is determined by" \
  " the \"-maxval\", \"-badOut\" and \"-isMaskOut\" arguments.  See below.\n" \
  "\n" \
  "  In any case, some filters interpret each input sample value {us} as a fuzzy" \
  " quantity, with a uniform distribution in an interval {[us-eps_us+eps]}, where" \
  " {eps} is determined from the input file's quantization step or" \
  " from a user-given \"-noise\" parameter.  For instance, the \"deviation\" filter" \
  " always returns a small positive value, even if all pixels" \
  " are equal.\n" \
  "\n" \
  "INVALID SAMPLE VALUES\n" \
  "  The program allows a specific sample value in the input image to be declared" \
  " as meaning `undefined value', so that any pixels with that value are excluded" \
  " from the filter computations.  See" \
  " the \"-badIn\" option below.  Pixels outside the image domain may or" \
  " may not be considered `undefined' (see \"OPTIONS\" below).  Similarly," \
  " a specific output sample value" \
  " may be used to encode `undefined' or `invalid' filter results.  See" \
  " the \"-badOut\" option below.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -maxval {OMAXVAL} \n" \
  "    Defines the nominal maximum sample value for the output" \
  " image.  Must be an integer between 1" \
  " and " stringify(PNM_MAX_MAXVAL) ".  If not specified," \
  " defaults to the input maximum sample value {IMAXVAL}.\n" \
  "\n" \
  "  -noise {NOISE} \n" \
  "    Defines the standard deviation of the quantization errors and" \
  " other noise that may affect the pixel values (after" \
  " conversion to [0_1] scale).  If not specified," \
  " defaults to the deviation of the quantization noise for uniform random variable" \
  " namely {sqrt(1/12)/maxval}.\n" \
  "\n" \
  "  -badIn {IBAD} \n" \
  "    If this option is present, and {IBAD} is in the range" \
  " {0..IMAXVAL}, any input sample with value {IBAD} is assumed" \
  " to be `undefined', and its weight is set to 0.  In that" \
  " case, the program implicitly subtracts 1 from any input" \
  " sample value in the range {IBAD+1..IMAXVAL}, so that the" \
  " nominal input range becomes {0..IMAXVAL-1}.  Thus, for" \
  " example, given \"-badIn 2\" and {IMAXVAL = 6}, input" \
  " sample values 0 through 6 are interpreted" \
  " as 0, 1, undefined, 2, 3, 4, and 5, respectively.\n" \
  "\n" \
  "  -badOut {OBAD} \n" \
  "    If this option is present, and {OBAD} is in" \
  " the range {0..OMAXVAL}, the output sample value {OBAD}" \
  " is reserved to represent `undefined' or `invalid' filter results.  In" \
  " that case, valid filter results will be scaled to the range {0..OMAXVAL-1}" \
  " instead of {0..OMAXVAL}, and then results in the range {OBAD..OMAXVAL-1} will" \
  " be incremented by 1. Thus, for example, given \"-badOut 2\" and" \
  " {OMAXVAL = 6}, the valid filter results will be scaled so as to range" \
  " from 0 through 5, and output values 0 through 6 will mean results" \
  " 0, 1, undefined, 2, 3, 4, and 5, respectively.\n" \
  "\n" \
  "    If \"-badOut\" is omitted but \"-badIn\" is given," \
  " and the input and output ranges are the same ({OMAXVAL=IMAXVAL})," \
  " then {OBAD} defaults to {IBAD} --- that is, the" \
  " same `undefined value' convention will be used for" \
  " the input and output images.  Otherwise," \
  " if \"-badOut\" is omitted, the output" \
  " will have no special `undefined' value." \
  " Valid filter results will be scaled to the" \
  " range {0..OMAXVAL}, and  invalid results will be mapped" \
  " to {OMAXVAL/2}, without any adjustment of higher values.\n" \
  "\n" \
  "  -keepBad \n" \
  "    This option forces an invalid output sample value {vs}" \
  " whenever the corresponding input sample {us} is" \
  " invalid.  It is effective" \
  " only if \"-badIn\" is specified, and only for some filters, such" \
  " as \"median\" and \"percentile\"," \
  " that might otherwise `repair' missing samples in the input.\n" \
  "\n" \
  "  -isMaskIn {IMK} \n" \
  "    This Boolean option specifies the method for conversion of input image" \
  " samples to float values in the real interval {[0_1]}.  Let {0..MAXVAL} be" \
  " the effective range of the integer samples from the input image, after" \
  " excluding the invalid value {IBAD} if given.  If the argument {IMK} is true (\"T\" or 1)," \
  " " sample_conv_0_1_isMask_true_INFO "  If the argument {IMK} is false (\"F\" or 0)," \
  " " sample_conv_0_1_isMask_false_INFO "  The default is \"-isMaskIn F\".\n" \
  "\n" \
  "  -isMaskOut {OMK} \n" \
  "    This Boolean option specifies the method for rounding computed float" \
  " samples, nominally in the real interval {[0_1]}, to integer samples" \
  " in the output image.  Let {0..MAXVAL} being the effective range of" \
  " the integer samples in the output image, namely {0..OMAXVAL} but" \
  " excluding the invalid value {OBAD} if given.  The two possibilities" \
  " selected by {OMK} are those described under the \"-isMaskIn\" option.  The" \
  " default is {OMK=IMK}.\n" \
  "\n" \
  "  -replicate \n" \
  "    This option specifies that any pixels of the" \
  " filter window that fall outside the image domain" \
  " should be assumed to have the same image value as the nearest pixel" \
  " in the domain.  In other words, the extremal rows" \
  " and columns are implicitly replicated outwards so" \
  " as to cover the window.  If not specified, pixels" \
  " that fall outsde the domain are assumed to be" \
  " undefined, and are given zero weight.\n" \
  "\n" \
  "  -excludeSelf \n" \
  "    This option forces the weight of the central window pixel {us}" \
  " to be zero, irrespective of its value in the window image file.\n" \
  "\n" \
  "  -filter {FKIND} {FPARAMS}.. \n" \
  "    Specifies the type of filter and its parameters (if any).  The" \
  " valid alternatives for {FKIND} and {FPARAMS} are:\n" \
  "\n" \
  PROG_INFO_FILTERS "\n" \
  "\n" \
  "    All the filter descriptions above assume that input" \
  " sample values and the {IMAXVAL} parameter have been" \
  " adjusted to exclude any `undefined' or `invalid' values," \
  " as described under the \"-badIn\" option.  The same" \
  " applies to the output values, if \"-badOut\" is in effect." \
  "\n" \
  "  -weights {WFILE}\n" \
  "    Specifies the name {WFILE} of a PGM or PBM image" \
  " file which contains the weights of pixels in the" \
  " window.  Its width and height must be" \
  " odd positive integers." \
  "\n" \
  "  -mask {MFILE}\n" \
  "    Specifies the name {MFILE} of a PGM or PBM image" \
  " file which contains a weight mask for the input image.  I must have the" \
  " same size as the input image.  If not specified, the" \
  " program assumes a mask of all ones --- that is, the pixel weights" \
  " are defined by the window image only.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  pnmnlfilt(1), pnmconvol(1), pgmenhance(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created with the name \"pgmwfilter\" on 2006-11-18" \
  " by J. Stolfi (IC-UNICAMP), based on Jef Poskanzer Netpbm image file" \
  " formats and earlier PGM filter programs.\n" \
  "\n" \
  "  The original median filter code was created by J. Stolfi on 1996-11-21," \
  " as a separate program \"pgmmedf\", inspired on the Netpbm tool \"pnmnlfilt\"" \
  " by Graeme W. Gill (1993 or earlier).\n" \
  "\n" \
  "  The original local normalization filter code was created by J. Stolfi" \
  " on 1996-11-21, as the separate program \"pgmlocnorm\".\n" \
  "\n" \
  "  The original rank filter code was created by J. Stolfi on 1999-04-18," \
  " as the separate program \"pgmrankf\".\n" \
  "\n" \
  "  The local stretch filter code was created by J. Stolfi on 2006-11-18.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2006-11-18 Created by J. Stolfi.\n" \
  "  2010-08-02 Renamed \"pnmwfilter\" and extended to arbitrary PNM input.\n" \
  "  2010-08-14 Added the \"-mask\" option.\n" \
  "  2010-08-14 Added the \"-iMaskIn\" and \"-iMaskOut\" options.\n" \
  "  2010-08-17 Added the \"variance\" filter.\n" \
  "  2010-08-17 Modified variance and deviation filter so that max is 1.0 for any noise.\n" \
  "  2010-08-17 Fixed (?) variance proc so that {fnoise} is the noise deviation in 0-1 scale.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmwfilter_C_COPYRIGHT "\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_FILTERS \
  "      median\n" \
  "        This filter sets {vs} to the median value of all" \
  " the samples in the window, taking their weights into" \
  " account.  The result {vs} is undefined if the total" \
  " weight of the window pixels is zero.\n" \
  "\n" \
  "      percentile {LEVEL}\n" \
  "        This filter is a generalization of the median" \
  " filter: it sets {vs} to the sample value that is greater than a" \
  " specified fraction {LEVEL} (in weight) of the" \
  " window pixels.  The result" \
  " is undefined if the total weight of window pixels is" \
  " zero.  Thus, \"percentile 0.5\"" \
  " gives the same effect as the median filter; \"percentile 0.0\" selects the" \
  " smallest sample with nonzero weight in the window; \"percentile 1.0\" selects the largest" \
  " such sample; \"percentile 0.25\" selects a sample that is greater than 25% of" \
  " the window samples; and so on.\n" \
  "\n" \
  "        In mathematical morphology nomenclature," \
  " \"percentile 0.0\" yields a grayscale erosion of the" \
  " image by the structuring element that consists of the window pixels with nonzero weight," \
  " while \"percentile 1.0\" yields a grayscale dilation by the same structuring element." \
  "\n" \
  "      average\n" \
  "        Replaces each sample {us} by the window average, namely the weighted average" \
  " of the samples in the window centered at {us}." \
  "\n" \
  "      variance\n" \
  "        Replaces each sample {us} by the window variance, defined as the weighted" \
  " average of the squares of the differences between the window average and the" \
  " samples in the window.  The square of the \"-noise\" parameter is added to the result." \
  "\n" \
  "      deviation\n" \
  "        Replaces each pixel {us} by the window standard deviation, defined as" \
  " the square root of the variance filter." \
  "      normalize {RMAG}\n" \
  "        This filter sets {vs} to the difference between" \
  " {us} and the weighted mean of all window pixels, relative" \
  " to the standard deviation of those pixels.  The result {vs} is" \
  " undefined if {us} is undefined, or if the total weight" \
  " of the window pixels is zero.\n" \
  "\n" \
  "        More precisely, this filter computes the average {M} and the standard" \
  " deviation {D} of all samples in the window.  The output" \
  " sample {vs} is {R = (us - M)/D}, linearly scaled" \
  " and shifted so that {R = -RMAG} maps to {vs = 0}, and {R = +RMAG}" \
  " maps to {vs = OMAXVAL}.  Thus, for example, {vs=0} means that {us}" \
  " is at least {RMAG} deviations below the window mean, {vs=OMAXVAL/2}" \
  " means {us} is equal to" \
  " the window's mean, and {vs=OMAXVAL} means {us} is at" \
  " least {RMAG} deviations above the mean.\n" \
  "\n" \
  "      stretch {LOLEVEL} {HILEVEL}\n" \
  "        This filter applies to the input sample {us} a linear scaling" \
  " and shift that would take the {LOLEVEL} percentile of the window samples" \
  " to value {LOLEVEL*OMAXVAL}, and the {HILEVEL} percentile to" \
  " value {HILEVEL*OMAXVAL}.  The result {vs} is" \
  " undefined if {us} is undefined, or if the total weight" \
  " of the window pixels is zero.\n" \
  "\n" \
  "        In particular, \"stretch 0.0 1.0\" sets {vs=0} if {us} is the lowest" \
  " window sample in the window, {vs=OMAXVAL} if {us}" \
  " is the largest sample in the window, and {vs=OMAXVAL/2} if {us}" \
  " is halfway between the largest and smallest window samples.\n" \
  "\n" \
  "      rank {LOLEVEL} {HILEVEL}\n" \
  "        This filter replaces {us} by its relative rank among the" \
  " window pixels, linearly scaled and shifted so that the relative" \
  " ranks {LOLEVEL} and {HILEVEL} map to 0 and {OMAXVAL}," \
  " respectively.  The result {vs} is" \
  " undefined if {us} is undefined, or if the total weight" \
  " of the window pixels is zero.\n" \
  "\n" \
  "        In particular, \"rank 0.0 1.0\" yields" \
  " {vs=0} if {us} is the smallest sample in the window," \
  " {vs=OMAXVAL} if {us} is the largest sample, and" \
  " {vs=OMAXVAL/2} if {us} is the median of the" \
  " window samples.\n" \
  "\n" \
  "        Note that the \"rank\" filter differs from the" \
  " \"stretch\" filter in the mapping of intermediate sample" \
  " values.  In standard image processing jargon, the" \
  " \"rank\" filter is a local version of `image equalization'," \
  " while \"stretch\" is a local version of `image normalization'." \

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <limits.h>
#include <values.h>

#include <jspnm.h> 
#include <jsfile.h> 
#include <jswsize.h> 
#include <ix.h> 
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <float_pnm_read_stream.h>
#include <float_pnm_write_stream.h>
#include <float_image_buffer.h>
#include <sample_conv.h>
#include <affirm.h> 
#include <argparser.h> 

#define INF INFINITY

/* DATA TYPES */

typedef enum
  { FT_KIND_AVERAGE,      /* Weighted average of the window samples */
    FT_KIND_VARIANCE,     /* Weighted variance of the window samples. */
    FT_KIND_DEVIATION,    /* Square root of variance. */
    FT_KIND_PERCENTILE,   /* Fixed percentile of window samples. */
    FT_KIND_RANK,         /* Rank of center sample among window samples. */
    FT_KIND_NORMALIZE,    /* Relative deviation of center sample from window mean. */
    FT_KIND_STRETCH       /* Relative position of center pixel between window percentiles. */
  } ft_kind_t;
  /* Kind of filter. */

typedef struct options_t 
  { char *iname;         /* Input filename ("-" for stdin). */
    char *mname;         /* Mask image name ("-" for stdin,{NULL} if not given). */ 
    char *wname;         /* Window/weight file name ("-" for stdin). */ 
    char *oname;         /* Output filename ("-" for stdout). */
    bool_t verbose;      /* TRUE says to mumble while working. */
    /* Input/output encoding and undefined value handling: */
    uint16_t maxval; /* Output maxval (0 if not specified). */
    uint32_t badIn;      /* Input sample value that means `undefined' ({>maxval} if none). */
    bool_t keepBad;      /* TRUE forces undef output when input is undef. */
    bool_t isMaskIn;     /* TRUE if input 0 and {maxval} are to be mapped to 0.0 and 1.0. */
    uint32_t badOut;     /* Output sample value that means `undefined' ({>maxval} if none). */
    bool_t isMaskOut;    /* TRUE if output 0 and {maxval} are to mean 0.0 and 1,0. */
    bool_t replicate;    /* TRUE specifies edge replication. */
    bool_t excludeSelf;  /* TRUE forces weight zero for center pixel. */
    /* Filter parameters (may be {NAN} if not applicable): */
    ft_kind_t ft_kind;   /* Kind of filter to apply. */
    double noise;        /* Standard deviation of assumed noise, or {NAN} if not given. */
    double level;        /* Level for percentile filter (0.5 for median). */
    double rmag;         /* Scaling coefficient for various filters. */
    double lolev, hilev; /* Low and high levels for stretch and rank filters. */
  } options_t;
  /* Arguments parsed from the command line. */

#define MAX_CHNS 3
  /* Max channels in a PNM file. */

typedef struct wpixel_t 
  { double val[MAX_CHNS]; /* Pixel samples from input image, scaled to {[0_1]}. */
    double msk;    /* Sample from mask image, scaled to {[0_1]}. */
    int32_t dx, dy;    /* Index displacements from window's refrence pixel. */
    double wht;    /* Weight from window image. */
  } wpixel_t;
  /* Data for a valid pixel extracted from the window, namely its
    samples {val[0..chns-1]}, the corresponding mask sample {msk},
    its position {dx,dy} relative to the window's
    reference pixel, and its window weight {wht}. 
    
    The fields {dx,dy,wht} are constant during processing,
    whereas {val} and {.msk} change as the window is scanned over the
    image. If {val[kc]} is {NAN} (meaning `undefined'), the pixel
    should be assumed to have zero weight, irrespective of the {wht}
    field. */

#define MAX_TOT_WEIGHT (((uint64_t)2)<<(50-16))
  /* Maximum sum of all raw (integer) weights in a window. Must be a
    power of 2, less than the precision of a {double}, less 16 bits because
    there may be mask weights also. */

#define MAX_WINDOW_SIZE 32767
  /* Maximum width and height of the window image. 
    Must be such that its square fits in an {int32_t}. */

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char* argv[]);

options_t *parse_options(int32_t argc, char **argv);

void read_window_image
  ( char *wname,         /* Window image's filename. */
    bool_t notSelf,      /* TRUE to exclude center pixel from window. */
    bool_t verbose,      /* TRUE to print diagnostics and statistics. */
    int32_t *wcolsP,         /* OUT: Effective window width (pixels). */
    int32_t *wrowsP,         /* OUT: Effective window height (pixels). */
    int32_t *nwP,            /* OUT: number of valid samples in window. */
    wpixel_t **wpixP     /* OUT: valid pixel values and weights in window. */
  );
  /* Reads the PGM or PBM image file called "{wname}" and creates the 
    corresponding window structure {*wpixP}.
    
    If {notSelf} is TRUE, pretends that the value of the center pixel is 0.
    
    In any case, collects all nonzero pixels into the {.wht} fields of
    an array of {wpixel_t} records, whose address is returned in
    {*wpixP} and whose length is returned in {*nwP}.
    
    The number of columns {wcols} (returned in {*wcolsP}) and the
    number of rows {wrows} (returned in {*wrowsP}) of the window image
    must be odd, and the central pixel in column {wcols/2} and
    {wrows/2} is assumed to be the window's reference pixel (with {.dx =
    .dy = 0}. The weights {.wht} are converted to {double} but
    otherwise not changed. If {verbose} is TRUE, prints various
    statistics about the window weights. */

void show_window_weights_statistics(FILE *wr, int32_t nw, wpixel_t wpix[], uint16_t maxval);
  /* Computes some statistics of the window weights as collected in
    {wpix[0..nw-1]} and writes them legibly to file {wr}.
    The {maxval} is used to convert weights to a {[0_1]} scale for
    redability. */

void filter_image_file
  ( FILE *ifile,
    FILE *mfile,
    FILE *ofile,
    options_t *o,
    int32_t wcols,
    int32_t wrows,
    int32_t nw,
    wpixel_t wpix[]
  );
  /* Reads a PBM, PGM, or PPM image from file {ifile}, applies the filter described
    by {o} with the window weights {wpix[0..nw-1].wht}, and writes the result to {ofile}.
    
    If {mfile} is not NULL, it must be a file containing a PGM or PBM mask image,
    with the same size as {ifile}'s, whose samples are intepreted as 
    weights of the corresponding {ifile} pixels.
    
    Assumes that the effective width of the window is {wcols} pixels
    and its effective height is {wrows}. Both must be odd. The
    window's reference pixel is assumed to be {wcols/2, wrows/2}. */

void grab_input_pixels_and_mask
  ( int32_t x,                      /* Column of current pixel in input image. */
    int32_t y,                      /* Row of current pixel in input image. */
    float_image_buffer_t *ibuf, /* Circular row buffer of input image. */
    float_image_buffer_t *mbuf, /* Circular row buffer of mask image, or NULL. */
    bool_t replicate,           /* TRUE extend image by edge replication. */
    int32_t nw,                     /* Number of valid samples in window. */
    wpixel_t wpix[]             /* IN/OUT: valid pixel values and weights in window. */
  );
  /* Updates the value field {val[0..chns-1]} of entries {wpix[0..nw-1]}, with
    data taken from channel {kc} of the input image, within the window
    centered at {x,y}. If {replicate} is FALSE, window pixels that
    fall outside the image's domain are set to {NAN}; otherwise they
    are copied from the nearest edge pixel.
    
    Assumes that the rows of the input image that are spanned by the
    window are stored in {buf}. The procedure also assumes that the fields
    {dx,dy,wht} of {wpix[0..nw-1]} have been set (see
    {read_window_image}), and it does not modify them. */

double apply_average_filter(int32_t kc, int32_t nw, wpixel_t wpix[]);
double apply_variance_filter(int32_t kc, int32_t nw, wpixel_t wpix[], double fnoise);
double apply_deviation_filter(int32_t kc, int32_t nw, wpixel_t wpix[], double fnoise);
double apply_percentile_filter(int32_t kc, int32_t nw, wpixel_t wpix[], int32_t wperm[], double fnoise, double level);
double apply_rank_filter(double ival, int32_t kc, int32_t nw, wpixel_t wpix[], double fnoise);
double apply_normalize_filter(double ival, int32_t kc, int32_t nw, wpixel_t wpix[], double fnoise, double rmag);
double apply_stretch_filter(double ival, int32_t kc, int32_t nw, wpixel_t wpix[], int32_t wperm[], double fnoise, double lolev, double hilev);
  /* These procedures apply the relevant filter to the channel {kc} samples 
    of the input image within the window centered at the current pixel, and
    return the output sample value for that pixel.
    
    These procedures assume that the channel samples, mask sample, and
    window weights of the pixels within the current window rectangle
    are stored in {wpix[0..nw-1].{val[kc],msk,wht}}. The effective
    weight of each pixel is the product {msk*wht}. 
    
    These procedures assume that no {val[kc]} is infinite, and they
    ignore pixels with {val[kc]==NAN} as if they had {msk==wht==0}.
    They assume that {msk} is never infinite or {NAN}. If the total
    weight of the valid pixels is zero, these procedures return {NAN}.
    
    The parameter {ival}, when required, should be the value of the
    current pixel in the input image. Note that this center pixel may
    or may not be present in {wpix}. */

void compute_weight_below_above(double ival, int32_t kc, int32_t nw, wpixel_t wpix[], double rv, double *sloP, double *shiP);
  /* Computes the sums {*sloP} and {*shiP} of the effective weight
    {.wht*.msk} for all valid pixels {wpix[0..nw-1]} that have values
    {.val[kc]} lower than higher, respectively, than {ival}. Actually
    assumes that the weight of each pixel {k} is spread over the
    interval of values {wpix[k].val[kc] ± rv}, so pixels that have
    values close enough to {ival} are split betwen the two sums. */

double compute_total_weight(int32_t kc, int32_t nw, wpixel_t wpix[]);
  /* Computes the sum of the effective weights {.wht*.msk} of all pixels {wpix[0..nw-1]}
    in the current window, ignoring any pixels with {.val[kc] == NAN}. */
    
double compute_window_avg(int32_t kc, int32_t nw, wpixel_t wpix[]);
  /* Computes the weighted mean of the sample values {wpix[0..nw-1].val[kc]}
    with the effective weights {.wht*.msk}.  Ignores pixels
    with {.val[kc] == NAN}; returns {NAN} if result is undefined. */
    
double compute_window_var(int32_t kc, int32_t nw, wpixel_t wpix[], double avg, double fnoise);
  /* Estimates the weighted mean squared deviation of the sample
    values {wpix[0..nw-1].val[kc]} from their weighted mean, with the
    effective weights {.wht*.msk}. Assumes that {wpix[0..nw-1].val[kc]}
    are the observed values, whose weighted average is {avg}. Assumes
    that the actual value of each sample is the given value plus a
    quantization noise with zero mean and deviation {fnoise}. Ignores
    pixels with {.val[kc] == NAN}. Always returns zero if there is only
    one pixel in the input list. Returns {NAN} if result is
    undefined. */

double find_percentile(int32_t kc, int32_t nw, wpixel_t wpix[], int32_t wperm[], double rv, double level);
  /* Finds the pixel value {fval} such that the total effective weight {.wht*.msk} of all
    window pixels {wpix[0..nv-1]} whose value {.val[kc]} is less than
    {fval} is {level} times sum {wtot} of all effective  weights.
    If the sum {wtot} is zero, returns {fval == level}.
    
    The permutation vector {wperm[0..nw-1]} must contains 
    the pixel indices {0..nw-1} sorted by increasing {val[kc]} field.
    Pixels with {.val[kc] == NAN} should be at the end of this list.
    
    Assumes that the true value of each valid pixel {wpix[k]} is
    actually {wpix[k].val[kc] + e} where {e} is an unknown error
    uniformly distributed in {-rv_+rv]}. So if that pixel is close
    enough to {fval}, it is counted with a fraction of its weight. */
    
int32_t compare_val(double a, double wa, double b, double wb);
  /* Compares two sample values in the order required by
    {wpixel_sort}. Namely, if {wa!=0} and {wb!=0}, returns {-1,0,+1} depending on whether
    {a<b}, {a==b}, or {a>b}, respectively.  As a special case,
    considers that a {NAN} value is greater than any other value, except that
    it is equal to {NAN}.  Also if {wa==0} treats {a} as {NAN},
    and if {wb==0} treats {b} as {NAN}. */

void wpixel_sort(int32_t kc, int32_t nw, wpixel_t wpix[], int32_t wperm[]);
  /* Assumes that {wperm[0..nw-1]} is a permutation of {0..nw-1}.
    Rearranges it so that {wpix[wperm[0..nw-1]].val[kc]} is
    non-decreasing. Entries with {NAN} sample values are
    sorted at the end of the list. */

/* ROUTINES */

int32_t main(int32_t argc, char* argv[])
  { /* Command line arguments: */
    options_t *o = parse_options(argc, argv);
    
    /* Consistency check: */
    assert((int32_t)PNM_MAX_SAMPLE < INT32_MAX - 1);

    /* Read input header and get image attributes: */
    FILE *ifile = open_read(o->iname, o->verbose);
    
    /* Write output file header: */
    FILE *ofile = open_write(o->oname, o->verbose);
    
    /* Get the relevant window pixels from the window image: */
    int32_t wcols; /* Window columns. */
    int32_t wrows; /* Window rows. */
    int32_t nw;
    wpixel_t *wpix;
    read_window_image(o->wname, o->excludeSelf, o->verbose, &wcols, &wrows, &nw, &wpix);
    
    /* Read the input image's mask file, if any: */
    FILE *mfile = (o->mname == NULL ? NULL : open_read(o->mname, o->verbose));

    filter_image_file(ifile, mfile, ofile, o, wcols, wrows, nw, wpix);
    
    if (ifile != stdin) { fclose(ifile); }
    if ((mfile != NULL) && (mfile != stdin)) { fclose(mfile); }
    if ((ofile != stdout) && (ofile != stderr)) { fclose(ofile); }
    return 0;
  }

void filter_image_file
  ( FILE *ifile,
    FILE *mfile,
    FILE *ofile,
    options_t *o,
    int32_t wcols,
    int32_t wrows,
    int32_t nw,
    wpixel_t wpix[]
  )
  { 
    /* Read the input file header, create a buffered read stream {istr} for it: */
    float_pnm_stream_t *istr = float_pnm_read_stream_new(ifile, o->isMaskIn, o->badIn, wrows);
    int32_t rows = istr->rows;
    int32_t cols = istr->cols;
    int32_t chns = istr->chns;
    uint16_t imaxval = istr->maxval;  /* Input image's {maxval}. */
    assert(imaxval > 0);
    assert((chns > 0) && (chns <= MAX_CHNS));
    float_image_buffer_t *ibuf = istr->buf;

    /* Determine the actual input sample maxval {irange}, skipping {o->badIn}: */
    uint16_t irange = (o->badIn <= (uint32_t)imaxval ? imaxval - 1 : imaxval);

    /* Get deviation of quantization noise in scaled input pixels: */
    double fstep = 1.0/(o->isMaskIn ? irange : irange+1);
    double fnoise = (! isnan(o->noise) ? o->noise : sqrt(1.0/12.0)*fstep);
    
    /* If a mask was given, create a stream {mstr} for reading it too: */
    float_pnm_stream_t *mstr;    /* Mask image's stream, or NULL. */
    uint16_t mmaxval;         /* Mask image's {maxval} (1 if none). */
    if (mfile != NULL)
      { mstr = float_pnm_read_stream_new(mfile, TRUE, PNM_NO_BADVAL, wrows); 
        demand(mstr->chns == 1, "mask image is not moochromatic");
        demand(mstr->cols == cols, "mask image has wrong width");
        demand(mstr->rows == rows, "mask image has wrong height");
        mmaxval = mstr->maxval;
        assert(mmaxval > 0);
      }
    else
      { mstr = NULL;
        mmaxval = 1;
      }
    float_image_buffer_t *mbuf = (mstr == NULL ? NULL : mstr->buf);
   
    /* Choose the output maxval: */
    uint16_t omaxval = (o->maxval != 0 ? o->maxval : imaxval);
    assert(omaxval > 0);
    
    /* Choose the output file parameters and write the output file header: */
    bool_t forceplain = FALSE;  /* TRUE to force the `plain' (ascii) format variant. */
    pnm_format_t oformat;       /* Format of output PNM file. */
    bool_t oraw;                /* TRUE if `raw' format variant. */
    bool_t obits;               /* TRUE for PBM format (`raw' or `plain'). */
    pnm_choose_output_format(imaxval, chns, forceplain,  &oformat, &oraw, &obits);
    pnm_write_header(ofile, cols, rows, omaxval, oformat);

    /* Allocate the output buffer for a single-row of quantized samples: */
    uint16_t *osmp = uint16_image_alloc_pixel_row(cols, chns);
    
    /* Check window dimensions: */
    demand(wcols % 2 == 1, "window width must be odd");
    demand(wrows % 2 == 1, "window height must be odd");
    int32_t wrady = wrows/2;  /* Half-height of window excluding center pixel. */
    
    /* Allocate and initialize the window permutation vectors: */
    int32_t *wperm[chns];
    for (int32_t kc = 0; kc < chns; kc++) 
      { wperm[kc] = notnull(malloc(nw*sizeof(int32_t)), "no mem");
        for (int32_t k = 0; k < nw; k++) { wperm[kc][k] = k; }
      }
   
    /* Loop on output image rows: */
    for (int32_t y = 0; y < rows; y++)
      { int32_t yimin = y - wrady; /* Min input {y} needed to compute this row. */
        if (yimin < 0) { yimin = 0; }
        int32_t yimax = y + wrady; /* Max input {y} needed to compute this row. */
        if (yimax >= rows) { yimax = rows-1; }
        
        /* Make sure that rows {yimin..yimax} are in the input buffer: */
        (void)float_pnm_read_stream_get_row(ifile, istr, yimax);
        assert(ibuf->yini <= yimin);
        assert(ibuf->ylim > yimax);
        
        /* Ditto for the mask buffer, if any: */
        if (mstr != NULL)
          { (void)float_pnm_read_stream_get_row(mfile, mstr, yimax);
            assert(mbuf->yini <= yimin);
            assert(mbuf->ylim > yimax);
          }
        
        /* Get samples {ismp[0..chna*cols-1]} of row {y} of input image: */
        double *ismp = float_image_buffer_get_row(ibuf, y);
        double *msmp = (mstr == NULL ? NULL : float_image_buffer_get_row(mbuf, y));
        
        /* Compute row {y} of output image: */
        for (int32_t x = 0; x < cols; x++)
          { /* Get the mask value for this input pixel: */
            double mval = (msmp == NULL ? 1 : msmp[x]);
            /* Get input pixel values and mask (if any) for the new window placement: */
            grab_input_pixels_and_mask(x, y, ibuf, mbuf, o->replicate, nw, wpix);
            /* Now loop on channels: */
            for (int32_t kc = 0; kc < chns; kc++)
              { /* Get the center sample {ival} of the input image (may be undefined): */
                double ival = (mval == 0 ? NAN : ismp[x*chns + kc]);
                /* Compute the floated output pixel value: */
                double oval = NAN;
                switch(o->ft_kind)
                  { case FT_KIND_AVERAGE:
                      oval = apply_average_filter(kc, nw, wpix);
                      break;
                    case FT_KIND_VARIANCE:
                      oval = apply_variance_filter(kc, nw, wpix, fnoise);
                      break;
                    case FT_KIND_DEVIATION:
                      oval = apply_deviation_filter(kc, nw, wpix, fnoise);
                      break;
                    case FT_KIND_PERCENTILE:
                      oval = apply_percentile_filter(kc, nw, wpix, wperm[kc], fnoise, o->level);
                      break;
                    case FT_KIND_RANK:
                      oval = apply_rank_filter(ival, kc, nw, wpix, fnoise);
                      break;
                    case FT_KIND_NORMALIZE:
                      oval = apply_normalize_filter(ival, kc, nw, wpix, fnoise, o->rmag);
                      break;
                    case FT_KIND_STRETCH:
                      oval = apply_stretch_filter(ival, kc, nw, wpix, wperm[kc], fnoise, o->lolev, o->hilev);
                      break;
                  }
                /* Quantize and store the output pixel value: */
                uint16_t oqts = pnm_quantize(oval, omaxval, o->isMaskOut, o->badOut);
                assert(oqts <= omaxval);
                osmp[x*chns + kc] = oqts;
              }
          }
        /* Write row {y} of the output image: */
        pnm_write_pixels(ofile, osmp, cols, chns, omaxval, oraw, obits);
      } 

    /* Release working storage: */
    for (int32_t kc = 0; kc < chns; kc++) { free(wperm[kc]); }
    float_pnm_stream_free(istr);
    if (mstr != NULL) { float_pnm_stream_free(mstr); }
    free(osmp);
  }

void read_window_image
  ( char *wname,         /* Window image's filename. */
    bool_t notSelf,      /* TRUE to exclude center pixel from window. */
    bool_t verbose,      /* TRUE to print diagnostics and statistics. */
    int32_t *wcolsP,     /* OUT: Effective window width (pixels). */
    int32_t *wrowsP,     /* OUT: Effective window height (pixels). */
    int32_t *nwP,        /* OUT: number of valid samples in window. */
    wpixel_t **wpixP     /* OUT: valid pixel values and weights in window. */
  )
  { 
    uint16_image_t *pim = uint16_image_read_pnm_named(wname, verbose);
    
    bool_t yup = FALSE; /* For now. */
    
    /* Get image dimensions: */
    if (pim->chns != 1) { pnm_error("window image file must be monochromatic"); }
    int32_t cols = pim->cols; 
    if (cols > MAX_WINDOW_SIZE) { pnm_error("window image too wide"); }
    if (cols % 2 != 1) { pnm_error("window image width must be odd"); }
    int32_t rows = pim->rows; 
    if (rows > MAX_WINDOW_SIZE) { pnm_error("window image too tall"); }
    if (rows % 2 != 1) { pnm_error("window image height must be odd"); }
    
    /* Effective half-width and half-height of non-zero pixels: */
    int32_t rx = 0;
    int32_t ry = 0;
    
    /* Max sample value in integer image: */
    uint16_t maxval = pim->maxval;
    
    /* Get the reference pixel (before Y flipping): */
    int32_t ctrx = cols / 2;
    int32_t ctry = rows / 2;

    /* Allocate window pixel list: */
    wpixel_t *wpix = (wpixel_t*)notnull(malloc(cols*rows*sizeof(wpixel_t)), "no mem");
    
    /* Basic counts: */
    int32_t nw = 0;         /* Number of nonzero weight samples. */
    int32_t nz = 0;      /* Number of zero weight samples. */
    uint64_t stot = 0;  /* Sum of all weight samples. */
  
    /* Save pixel values in {wpix[k].wht}, get {nw,nz,stot}: */
    int32_t x, y;
    for(y = 0; y < rows; y++)
      { /* Get the PNM row {pnmy} from the logical row {y}: */
        int32_t pnmy = (yup ? rows - 1 - y : y);
        uint16_t *prow = pim->smp[pnmy];
        for(x = 0; x < cols; x++)
          { /* Convert int32_t sample {*prow} to float {v}, store, keep stats: */
            uint16_t s = (*prow);
            assert(s <= maxval);

            /* Exclude center pixel if so specified: */
            if (notSelf && (x == ctrx) && (y == ctry))
              { if (verbose) { fprintf(stderr, "excluding pixel in col %d row %d\n", ctrx, ctry); }
                s = 0;
              }

            int32_t dx = x-ctrx, dy = y-ctry; /* Indices relative to center. */

            if (s == 0)
              { /* Count and ignore: */
                nz++;
              }
            else
              { /* Save weight and displacements: */
                wpixel_t *wp = &(wpix[nw]);
                wp->dx = dx; wp->dy = dy;
                wp->wht = (double)s;
                for (int32_t kc = 0; kc < MAX_CHNS; kc++) { wp->val[kc] = NAN; } /* Just in case. */
                wp->msk = 1;   /* Just in case. */
                nw++;
                /* Update true extent: */
                if (abs(dx) > rx) { rx = abs(dx); }
                if (abs(dy) > ry) { ry = abs(dy); }
                /* Update statistics: */
                if (s > MAX_TOT_WEIGHT - stot)
                  { fprintf(stderr, ("total weight in window file exceeds %" uint64_u_fmt "\n"), MAX_TOT_WEIGHT);
                    assert(FALSE);
                  }
                stot += s;
              }

            prow++;
          }
      }
      
    /* The image is not needed anymore: */
    uint16_image_free(pim);
    
    demand(nw > 0, "window image is completely zero");

    /* Trim excess allocated space: */
    wpix = (wpixel_t*)notnull(realloc(wpix, nw*sizeof(wpixel_t)), "no mem");

    /* Compute effective window size and center: */
    int32_t wcols = 2*rx + 1;
    int32_t wrows = 2*ry + 1;
    int32_t wcx = wcols/2;
    int32_t wcy = wrows/2;
         
    if (verbose) 
      { /* Print statistics: */
        int32_t npixels = cols*rows;
        int32_t wnp = wcols*wrows;
        fprintf(stderr, "window image statistics\n");
        fprintf(stderr, "  raw window image dimensions %d × %d = %d pixels\n", cols, rows, npixels);
        fprintf(stderr, "  effective window dimensions %d × %d = %d pixels \n", wcols, wrows, wnp);
        fprintf(stderr, "  window rectangle [%+d..%+d] × [%+d..%+d]\n", 0-wcx, wcols-wcx, 0-wcy, wrows-wcy);
        fprintf(stderr, "  %10d zero samples\n", nz);
        fprintf(stderr, "  %10d nonzero samples\n", nw);
        fprintf(stderr, ("  %10" uint64_u_fmt " total raw weight\n"), stot);
        show_window_weights_statistics(stderr, nw, wpix, maxval);
      }

    /* Return results: */
    *wcolsP = wcols;
    *wrowsP = wrows;
    *nwP = nw;
    *wpixP = wpix;
  }

void show_window_weights_statistics(FILE *wr, int32_t nw, wpixel_t wpix[], uint16_t maxval)
  { 
    double scale = (double)maxval;
    double min_w = +INF;  /* Min nonzero weight. */
    double max_w = -INF;  /* Max nonzero weight. */
    int32_t nmin = 0;
    int32_t nmax = 0;
    double ext[2] = {0, 0}; /* Max extent in each axis. */
    double rad = 0;         /* Max radial extent. */
    double sum_w = 0;       /* Sum of weights. */
    double sum_w2 = 0;      /* Sum of weights squared. */
    double sum_wx = 0;      /* Sum of weights times X coordinate. */
    double sum_wy = 0;      /* Sum of weights times Y coordinate. */
    double sum_wx2 = 0;     /* Sum of weights times X^2. */
    double sum_wy2 = 0;     /* Sum of weights times Y^2. */
    for (int32_t k = 0; k < nw; k++)
      { wpixel_t *wp = &(wpix[k]);
        double w = wp->wht;
        double dx = wp->dx;
        double dy = wp->dy;
        if (w < min_w) { min_w = w; nmin = 0; }
        if (w == min_w) { nmin++; }
        if (w > max_w) { max_w = w; nmax = 0; }
        if (w == max_w) { nmax++; }
        sum_w += w;
        sum_w2 += w*w;
        sum_wx += w*dx;
        sum_wy += w*dy;
        /* Assume that weight is uniformly distributed in pixel: */
        sum_wx2 += w*(dx*dx + 1.0/12.0);
        sum_wy2 += w*(dy*dy + 1.0/12.0);
        /* Update sample extents and radius: */
        ext[0] = fmax(ext[0], fabs(dx) + 0.5);
        ext[1] = fmax(ext[1], fabs(dy) + 0.5);
        rad = fmax(rad, hypot(fabs(dx)+0.5, fabs(dy)+0.5));
      }
    if (sum_w == 0) { return; }
    
    /* Compute barycenter relative to nominal center: */
    double ctr[2];
    ctr[0] = sum_wx/sum_w;
    ctr[1] = sum_wy/sum_w;
            
    /* Compute sample moments around nominal center: */
    double mmt[2];
    mmt[0] = sqrt(sum_wx2/sum_w);
    mmt[1] = sqrt(sum_wy2/sum_w);
    double tot_mmt = sqrt((sum_wx2 + sum_wy2)/sum_w);
    
    fprintf(wr, "  min weight =   %12.0f = %17.15f (×%d)\n", min_w, min_w/scale, nmin);
    fprintf(wr, "  max weight =   %12.0f = %17.15f (×%d)\n", max_w, max_w/scale, nmax);
    fprintf(wr, "  tot weight =   %12.6f = %17.15f \n", sum_w, sum_w/scale);
    fprintf(wr, "  barycenter =   ( %8.3f %8.3f ) pixels\n", ctr[0], ctr[1]);
    fprintf(wr, "  max extents =  ( %8.3f %8.3f ) pixels\n", ext[0], ext[1]);
    fprintf(wr, "  max radius =   %10.3f pixels\n", rad);
    fprintf(wr, "  axis moments = ( %8.3f %8.3f ) pixels\n", mmt[0], mmt[1]);
    fprintf(wr, "  total moment = %8.3f pixels\n", tot_mmt);
  }

void grab_input_pixels_and_mask
  ( int32_t x,                      /* Column of current pixel in input image. */
    int32_t y,                      /* Row of current pixel in input image. */
    float_image_buffer_t *ibuf, /* Circular row buffer of input image. */
    float_image_buffer_t *mbuf, /* Circular row buffer of mask image. */
    bool_t replicate,           /* TRUE extend image by edge replication. */
    int32_t nw,                     /* Number of valid samples in window. */
    wpixel_t wpix[]             /* IN/OUT: valid pixel values and weights in window. */
  )
  { int32_t NC = ibuf->sz[0];
    int32_t NX = ibuf->sz[1];
    int32_t NY = ibuf->sz[2];
    if (mbuf != NULL)
      { assert(mbuf->sz[0] == 1);
        assert(mbuf->sz[1] == NX);
        assert(mbuf->sz[2] == NY);
      }
   
    /* Scan the window pixels with non-zero weight: */
    for (int32_t kw = 0; kw < nw; kw++)
      { /* Get the pixel's entry in {wpix} table: */
        wpixel_t *wp = &(wpix[kw]);
        /* Compute pixel's indices {xp,yp} in input img, or -1 if non-existant: */
        int32_t yp = y + wp->dy;
        int32_t xp = x + wp->dx;
        if (replicate)
          { /* Map invalid {xp,yp} to nearest pixel in domain: */
            if (yp < 0) { yp = 0; } else if (yp >= NY) { yp = NY - 1; }
            if (xp < 0) { xp = 0; } else if (xp >= NX) { xp = NX - 1; }
          }
        else
          { /* Set invalid {xp,yp} to -1: */
            if ((xp < 0) || (xp >= NX) || (yp < 0) || (yp >= NY)) { xp = yp = -1; }
          }
        /* Get and set the value: */  
        if (yp < 0)
          { /* Pixel doesn't exist: */
            for (int32_t kc = 0; kc < NC; kc++) { wp->val[kc] = NAN; }
            wp->msk = 0;
          }
        else
          { /* Pixel exists, get sample values from input image: */
            double *ismp = float_image_buffer_get_row(ibuf, yp);
            for (int32_t kc = 0; kc < NC; kc++) { wp->val[kc] = ismp[xp*NC + kc]; }
            if (mbuf == NULL)
              { wp->msk = 1; }
            else
              { double *msmp = float_image_buffer_get_row(mbuf, yp);
                wp->msk = msmp[xp];
              }
          }
      }
  }

/* FILTERS */

#define M_SQRT3 (1.73205080756887729352)
   /* Radius of range of a uniform random variable with unit deviation. */

double apply_average_filter(int32_t kc, int32_t nw, wpixel_t wpix[])
  { return compute_window_avg(kc, nw, wpix); }

double apply_variance_filter(int32_t kc, int32_t nw, wpixel_t wpix[], double fnoise)
  { double avg = compute_window_avg(kc, nw, wpix);
    double var = compute_window_var(kc, nw, wpix, avg, fnoise);
    /* Return deviation, scaled so that the maximum value is 1: */
    return var/(0.25 + fnoise*fnoise);
  }

double apply_deviation_filter(int32_t kc, int32_t nw, wpixel_t wpix[], double fnoise)
  { /* Compute variance, scaled so that the maximum value is 1: */
    double var = apply_variance_filter(kc, nw, wpix, fnoise);
    return sqrt(var);
  }
 
double apply_percentile_filter(int32_t kc, int32_t nw, wpixel_t wpix[], int32_t wperm[], double fnoise, double level)
  { /* Sort the pixels by value: */
    wpixel_sort(kc, nw, wpix, wperm);
    /* Find the value that is the percentile at fraction {level}: */
    double rv = M_SQRT3*fnoise; /* Width of each interval in value axis: */
    return find_percentile(kc, nw, wpix, wperm, rv, level);
  }

double apply_normalize_filter(double ival, int32_t kc, int32_t nw, wpixel_t wpix[], double fnoise, double rmag)
  { if (isnan(ival)) { return NAN; }
    /* Compute the weighted average and deviation of pixel values: */
    double avg = compute_window_avg(kc, nw, wpix);
    if (isnan(avg)) { return NAN; }
    double dev = sqrt(compute_window_var(kc, nw, wpix, avg, fnoise));
    if (isnan(dev)) { return NAN; }
    /* Rescale with {avg - rmag*dev} as black, {avg  + rmag*dev} as white: */
    if (dev == 0.0) { return NAN; }
    double rval = (ival - avg)/dev;
    return (rval/rmag + 1)/2;
  }
  
double apply_stretch_filter(double ival, int32_t kc, int32_t nw, wpixel_t wpix[], int32_t wperm[], double fnoise, double lolev, double hilev)
  { if (isnan(ival)) { return NAN; }
    /* Sort the pixels by value: */
    wpixel_sort(kc, nw, wpix, wperm);
    /* Find the value {loval} of the low percentile at level {lolev}: */
    double rv = M_SQRT3*fnoise; /* Width of each interval in value axis: */
    double loval = find_percentile(kc, nw, wpix, wperm, rv, lolev);
    /* Find the value {hival} of the high percentile at level {hilev}: */
    double hival = (lolev == hilev ? loval : find_percentile(kc, nw, wpix, wperm, rv, hilev));
    /* Find the middle point: */
    double mdlev = (lolev + hilev)/2;
    double mdval = (loval + hival)/2;
    /* Find the position of {ival} relative to {loval,hival}: */
    double scale;
    if (loval == hival)
      { if (lolev == hilev)
          { scale = 1.0; }
        else
          { return NAN; }
      }
    else
      { scale = (hilev - lolev)/(hival - loval); }
    return mdlev + scale*(ival - mdval);
  }

double apply_rank_filter(double ival, int32_t kc, int32_t nw, wpixel_t wpix[], double fnoise)
  {
    if (isnan(ival)) { return NAN; }
    /* Compute total weight of pixels above and below {ival} : */
    double rv = M_SQRT3*fnoise; /* Width of each interval in value axis: */
    double sum_lo, sum_hi;
    compute_weight_below_above(ival, kc, nw, wpix, rv, &sum_lo, &sum_hi);
    /* Compute pixel rank: */
    double tot = sum_lo + sum_hi;
    if (tot == 0.0)
      { return NAN; }
    else
      { return sum_lo/tot; }
  }

/* AUXILIARY PROCS */

void compute_weight_below_above(double ival, int32_t kc, int32_t nw, wpixel_t wpix[], double rv, double *sloP, double *shiP)
  {
    double sum_lo = 0; /* Weight below {ival}. */
    double sum_hi = 0; /* Weight above {ival}. */
    for (int32_t k = 0; k < nw; k++) 
      { wpixel_t *wp = &(wpix[k]);
        double s = wp->val[kc], w = wp->wht*wp->msk;
        if (! isnan(s))
          { double ssup = s + rv;
            double sinf = s - rv;
            if (ssup <= ival) 
              { sum_lo += w; }
            else if (sinf >= ival)
              { sum_hi += w; }
            else
              { /* Split fuzzy pixel betwen low and high: */
                double rlo = (ival - sinf)/(2*rv), rhi = 1 - rlo;
                sum_lo += rlo*w;
                sum_hi += rhi*w;
              }
          }
      }
    (*sloP) = sum_lo; (*shiP) = sum_hi;
  }

double compute_window_avg(int32_t kc, int32_t nw, wpixel_t wpix[])
  { if (nw == 0) { return NAN; }
    double sum_w = 0.0, sum_sw = 0.0;
    for (int32_t k = 0; k < nw; k++) 
      { wpixel_t *wp = &(wpix[k]);
        double s = wp->val[kc];
        if (! isnan(s))
          { double w = wp->wht*wp->msk/((double)MAX_TOT_WEIGHT);
            sum_w += w;
            sum_sw += s*w;
          }
      }
    return (sum_w == 0.0 ? NAN : sum_sw/sum_w);
  }

double compute_window_var(int32_t kc, int32_t nw, wpixel_t wpix[], double avg, double fnoise)
  { if (isnan(avg)) { return NAN; }
    if (nw == 0) { return NAN; }
    if (nw == 1) { return 0; }
    double sum_w = 0.0;
    double sum_w2 = 0;
    double sum_wd2 = 0.0;
    for (int32_t k = 0; k < nw; k++) 
      { wpixel_t *wp = &(wpix[k]);
        double s = wp->val[kc];
        if (! isnan(s))
          { double d = s - avg;
            double w = wp->wht*wp->msk/((double)MAX_TOT_WEIGHT);
            sum_w += w; 
            sum_w2 += w*w; 
            sum_wd2 += w*d*d;
          }
      }
    if (sum_w == 0.0) { return NAN; }
    /* Contribution of noise: */
    double fnoise2 = fnoise*fnoise;
    return sum_wd2/sum_w + fnoise2;
  }
  
double find_percentile(int32_t kc, int32_t nw, wpixel_t wpix[], int32_t wperm[], double rv, double level)
  { 
    double wtot = compute_total_weight(kc, nw, wpix);
    if (wtot == 0) { return level; }
    
    /* Assume that {wperm} is sorted so that {wpix[wperm[0..nw]].val[kc]} is
      non-increasing, with {NAN}s at the end. Reduce {nw} to
      exclude {NAN} pixels: */
    while ((nw > 0) && isnan(wpix[wperm[nw-1]].val[kc])) { nw--; }
    if (nw == 0) { return level; }
    assert(isfinite(wpix[wperm[0]].val[kc]));
    assert(isfinite(wpix[wperm[nw-1]].val[kc]));
    
    /* Useful parameters: */
    double dens = 1.0/(2*rv);   /* Weight density inside an interval of unit weight. */
    
    /* Trivial cases: */
    if (level <= 0.0) { return wpix[wperm[0]].val[kc] - rv; }
    if (level >= 1.0) { return wpix[wperm[nw-1]].val[kc] + rv; }

    /* Define search direction parameters: */
    int32_t dir;      /* Direction of scan. */
    int32_t ini, fin; /* First and last indices in search dir. */
    if (level <= 0.5)
      { dir = +1; ini = 0; fin = nw-1; }
    else
      { dir = -1; ini = nw-1; fin = 0; level = 1.0 - level; }
    assert((level >= 0) && (level <= 0.5));

    /* Compute the total weight below the desired percentile level: */
    double wgoal = level*wtot;
    
    /* First and last value of {v} in search order: */
    double vini = wpix[wperm[ini]].val[kc] - dir*rv;
    double vfin = wpix[wperm[fin]].val[kc] + dir*rv;
    
    /* Search state: */
    double v;            /* Candidate for percentile value. */
    double w;            /* Total weight between {vini} and {v}. */
    int32_t km, kp;          /* First and last intervals that overlap {v}, in search order. */
    double wm, wp;       /* Total weight of pixels up to {km-1} and up to {kp-1}. */
    int32_t kkm, kkp;        /* Ditto, actual indices. */
    wpixel_t *wkm, *wkp; /* Pointers to {wpix[wperm[kkm]],wpix[wperm[kkp]]}. */
    /* State invariant:
        (0) {kkp == ini + kp*dir} and {kkm == ini + km*dir}.
        (9) {wkp == &(wpix[wperm[kkp]])} and {wkm == &(wpix[wperm[kkm]])}.
        (1) {v == wpix[wperm[k]].val ± rv} for some {k} in {0..nw-1}.
        (2) {km} is the min index in {0..nw-1} such that {(v - wkm->val[kc])*dir <= rv}.
        (3) {kp} is the min index in {0..nw-1} such that {(wkp->val[kc] - v)*dir > rv}, or {nw}.
        (4) {wm} is the total weight of pixels from {ini} to {kkm} including {kkm}.
        (5) {wp} is the total weight of pixels from {ini} to {kkp} excluding {kkp}.
        (6) {w} is the total weight up to {v} (including partial intervals).
      therefore the interval {wpix[wperm[ini + k*dir]].val[kc] ± rv} contains {v} iff {k \in {km..kp-1}}.
      Note that always {km < kp}.
    */
    /* Seach for {v} in direction {dir} until {w >= wgoal}: */
    v = vini; w = 0;
    km = 0; kkm = ini; wm = 0; wkm = &(wpix[wperm[kkm]]);
    kp = 0; kkp = ini; wp = 0; wkp = &(wpix[wperm[kkp]]);
    do
      { assert(isfinite(wkp->val[kc]));
        wp += wkp->wht*wkp->msk; 
        kp++; kkp += dir; wkp = &(wpix[wperm[kkp]]);
      }
    while((kp < nw) && (wkp->val[kc] - dir*rv == v));

    while ((w < wgoal) && (v != vfin))
      { /* Paranoia checks: */
        assert((0 <= km) && (km < nw));
        assert((0 < kp) && (kp <= nw));
        assert(km < kp);
        /* Exit all intervals that end at {v}: */
        while (TRUE)
          { if (v != wkm->val[kc] + dir*rv) { break; }
            wm += wkm->wht*wkm->msk; 
            km++; kkm += dir; wkm = &(wpix[wperm[kkm]]); 
            assert(km < nw);
          }
        /* Get the two candidates for the next {v}: */
        double evm = wkm->val[kc] + dir*rv;
        assert((evm - v)*dir > 0); 
        double evp = (kp >= nw ? ((double)dir)*INF : (double)(wkp->val[kc] - dir*rv));
        if (! ((evp - v)*dir > 0)) 
          { fprintf(stderr, "evp = %23.15e  v = %23.15e dir = %+2d\n", evp, v, dir);
            for (int32_t i = 0; i < nw; i++)
              { wpixel_t *wi = &(wpix[wperm[i]]); 
                fprintf(stderr, "  %04d  dx = %+3d  dy = %+3d", i, wi->dx, wi->dy);
                fprintf(stderr, "  val = %23.15e msk = %23.15e wht = %23.15e\n", wi->val[kc], wi->msk, wi->wht);
              }
          }
        assert((evp - v)*dir > 0); 
        /* Choose the earliest one: */
        double vnext = ((evm - evp)*dir <= 0 ? evm : evp);
        /* Compute the total weight of the intervals that cover the gap to {vnext}: */
        double wspan = wp - wm;
        /* Compute the total weight {wnext} up to {vnext}: */
        double wnext = wm + wspan*dens*fabs(vnext - v);
        /* Make sure that {w} always increases, in spite of roundoff: */
        if (wnext < w) { wnext = w; }
        /* Are we done? */
        if (wnext < wgoal)
          { /* The percentile is beyond {vnext}, move {v} to {vnext}: */
            v = vnext; w = wnext;
            /* Enter all intervals that begin at {v}: */
            while((kp < nw) && (v == wkp->val[kc] - dir*rv))
              { assert(isfinite(wkp->val[kc]));
                wp += wkp->wht*wkp->msk; 
                kp++; kkp += dir; wkp = &(wpix[wperm[kkp]]);
              }
          }
        else
          { /* The percentile lies between {v} and {vnext}, interpolate: */
            double dw = wnext - w;
            double hw = wgoal - w;
            if (dw > 0) { v = v + (vnext - v)*(hw/dw); }
            w = wgoal; /* This should force the iteration to stop. */
          }
      }
      
    /* Here we should have {w == wgoal} except perhaps for roundoff losses. */
    /* Hence {v} is the desired percentile: */
    return v;
  }
 
double compute_total_weight(int32_t kc, int32_t nw, wpixel_t wpix[])
  { double wtot = 0; 
    for (int32_t k = 0; k < nw; k++) 
      { wpixel_t *wkp = &(wpix[k]); 
        if (! isnan(wkp->val[kc]))
          { double wk = wkp->wht*wkp->msk;
            double wnext = wtot + wk;
            /* The weights must be summable without roundoff error: */
            assert(wnext - wtot == wk); 
            wtot = wnext;
          }
      }
    return wtot;
  }

options_t *parse_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
    
    /* Get keyword arguments: */

    o->verbose = argparser_keyword_present(pp, "-verbose");
    
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (int16_t)argparser_get_next_int(pp, 0, PNM_FILE_MAX_MAXVAL); }
    else
      { /* Use input maxval: */
        o->maxval = 0;
      }
 
    if (argparser_keyword_present(pp, "-badIn"))
      { o->badIn = (uint32_t)argparser_get_next_int(pp, 0, UINT32_MAX); }
    else
      { /* No undef input value: */
        o->badIn = PNM_NO_BADVAL;
      }
 
    o->keepBad = argparser_keyword_present(pp, "-keepBad");
    
    if (argparser_keyword_present(pp, "-badOut"))
      { o->badOut = (uint32_t)argparser_get_next_int(pp, 0, UINT32_MAX); }
    else
      { /* No undef output value: */
        o->badOut = PNM_NO_BADVAL;
      }
 
    if (argparser_keyword_present(pp, "-isMaskIn"))
      { o->isMaskIn = argparser_get_next_bool(pp); }
    else
      { /* No undef input value: */
        o->isMaskIn = FALSE;
      }
 
    if (argparser_keyword_present(pp, "-isMaskOut"))
      { o->isMaskOut = argparser_get_next_bool(pp); }
    else
      { /* No undef input value: */
        o->isMaskOut = o->isMaskIn;
      }
 
    o->replicate = argparser_keyword_present(pp, "-replicate");
    
    o->excludeSelf = argparser_keyword_present(pp, "-excludeSelf");
    
    /* Quantization noise kind and magnitude: */
    if (argparser_keyword_present(pp, "-noise"))
      { o->noise = argparser_get_next_double(pp, 0.0, 1.0e+5); }
    else
      { /* Compute it later: */
        o->noise = NAN;
      } 

    /* Window image: */
    argparser_get_keyword(pp, "-weights");
    o->wname = argparser_get_next(pp);
    
    /* Mask image: */
    if (argparser_keyword_present(pp, "-mask"))
      { o->mname = argparser_get_next(pp); }
    else
      { /* Standard deviation of a uniform variable in {[-0.5 _ +0.5]}: */
        o->mname = NULL;
      } 

    /* Filter kind and parameters: */
    o->level = NAN;
    o->lolev = NAN;
    o->hilev = NAN;
    o->rmag = NAN;
    if (argparser_keyword_present(pp, "-filter"))
      { 
        if (argparser_keyword_present_next(pp, "average"))
          { o->ft_kind = FT_KIND_AVERAGE; }
        else if (argparser_keyword_present_next(pp, "deviation"))
          { o->ft_kind = FT_KIND_DEVIATION; }
        else if (argparser_keyword_present_next(pp, "variance"))
          { o->ft_kind = FT_KIND_VARIANCE; }
        else if (argparser_keyword_present_next(pp, "percentile"))
          { o->ft_kind = FT_KIND_PERCENTILE;
            o->level = argparser_get_next_double(pp, 0.0, 1.0);
          }
        else if (argparser_keyword_present_next(pp, "median"))
          { o->ft_kind = FT_KIND_PERCENTILE;
            o->level = 0.5;
          }
        else if (argparser_keyword_present_next(pp, "rank"))
          { o->ft_kind = FT_KIND_RANK; }
        else if (argparser_keyword_present_next(pp, "normalize"))
          { o->ft_kind = FT_KIND_NORMALIZE;
            o->rmag = argparser_get_next_double(pp, 1.0e-6, 1.0e+6);
          }
        else if (argparser_keyword_present_next(pp, "stretch"))
          { o->ft_kind = FT_KIND_STRETCH;
            o->lolev = argparser_get_next_double(pp, -1.0e+6, +1.0e+6);
            o->hilev = argparser_get_next_double(pp, -1.0e+6, +1.0e+6);
          }
        else 
          { argparser_error(pp, "invalid filter kind"); }
      }
    else
      { argparser_error(pp, "missing \"-filter\" option"); }

    /* Skip to positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Get optional input file name: */
    if (argparser_next(pp) != NULL)
      { o->iname = argparser_get_next(pp); }
    else
      { o->iname = "-"; }
 
    /* Get optional output file name: */
    if (argparser_next(pp) != NULL)
      { o->oname = argparser_get_next(pp); }
    else
      { o->oname = "-"; }
    
    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }

int32_t compare_val(double a, double wa, double b, double wb)
{ if ((wa == 0) || isnan(a))
    { return ((wb == 0) || isnan(b) ? 0 : +1); }
    else if ((wb == 0) || isnan(b))
      { return -1; }
    else if (a < b)
      { return -1; }
    else if (a > b)
      { return +1; }
    else
      { return 0; }
  }
  
void wpixel_sort(int32_t kc, int32_t nw, wpixel_t wpix[], int32_t wperm[])
  { 
    auto int32_t compare_indices(const void *akp, const void *bkp);
    /* Given two pointers to two elements of {wperm[0..nw-1]}, compares
       the corresponding values of {wpix[].val[kc]}  with {compare_val}. */
    
    int32_t compare_indices(const void *akp, const void *bkp)
      { int32_t ak = *((int32_t *)akp); assert((ak >= 0) && (ak < nw)); wpixel_t *wap = &(wpix[ak]);
        int32_t bk = *((int32_t *)bkp); assert((bk >= 0) && (bk < nw)); wpixel_t *wbp = &(wpix[bk]);
        return compare_val(wap->val[kc], wap->wht*wap->msk, wbp->val[kc], wbp->wht*wbp->msk);
      }
  
    qsort(wperm, nw, sizeof(int32_t), compare_indices);
    /* Paranoia check: */
    for (int32_t i = 1; i < nw; i++) 
      { wpixel_t *wap = &(wpix[wperm[i-1]]);
        wpixel_t *wbp = &(wpix[wperm[i]]);
        assert(compare_val(wap->val[kc], wap->wht*wap->msk, wbp->val[kc], wbp->wht*wbp->msk) <= 0);
      }
  }
