#define PROG_NAME "pnmdiffop"
#define PROG_DESC "PBM/PGM/PPM differential operatios -- gradient,hessian,laplacian,elongation"
#define PROG_VERS "1.0"

#define pnmdiffop_C_COPYRIGHT \
  "Copyright © 2006 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-09-24 12:00:47 by stolfi */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -op { " image_window_op_HELP_OPERATORS " } \\\n" \
  "  [ -smoothed ] \\\n" \
  "  [ -squared ] \\\n" \
  "  [ -scale {SCALE} ] \\\n" \
  "  [ -offset {OFFSET} ] \\\n" \
  "  [ -maxval {OMAXVAL} ] \\\n" \
  "  [ -badIn {IBAD} [ -keepBad ] ] \\\n" \
  "  [ -badOut {OBAD} ] \\\n" \
  "  [ -isMaskIn {IMK} ] \\\n" \
  "  [ -isMaskOut {OMK} ] \\\n" \
  "  [ -replicate ] \\\n" \
  "  [ -yup ] \\\n" \
  "  [<] INFILE \\\n" \
  "  [>] OUTFILE "

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program applies one of several local operator to an image {U}, read" \
  " from the PBM, PGM or PPM file {INFILE}.  The operator is" \
  " selected by the \"-op\", \"-smoothed\", and \"-square\" options.  The output is another" \
  " image {V} with the same number of channels and same dimensions, written to the" \
  " file {OUTFILE} with the same sample" \
  " type (PBM, PGM, or PPM) and same encoding (`raw' or `plain').\n" \
  "\n" \
  "  For monochromatic" \
  " images, the sample of {V} at any pixel is" \
  " computed from the samples of {U} within a 3 by 3 window centered at the" \
  " same pixel of {U}.  For color (PPM) images," \
  " each channel is processed independently, as a separate monochromatic image.\n" \
  "\n" \
  "IMAGE ENCODING\n" \
  "  The samples in the input and output images are assumed to have linear" \
  " encoding (gamma 1.0).\n" \
  "\n" \
  "   First, the input image samples are converted to float values" \
  " in {[0_1]}, as determined by the \"-badIn\" and \"-isMaskIn\" arguments.\n" \
  "\n" \
  "  Then the selected operator is applied to the 3 by 3 window centered at each" \
  " pixel, giving a sample value for the corresponding pixel of the output image.  This" \
  " value will lie in some /natural range/ characteristic of the operator, that" \
  " is either {[0 _ r]} or {[-r _ +r]} for some positive value {r}.  See the" \
  " section OPERATOR RANGES below.\n" \
  "\n" \
  "  The operator result is then optionally multiplied by a given {SCALE} factor. See" \
  " the \"-scale\" option below.  Then the value is mapped by an affine function" \
  " that takes the natural range (ignoring the scaling) to the" \
  " interval {[0 _ 1]}.  Then an optional {OFFSET} is added to it; see" \
  " the \"-offset\" option below.\n" \
  "\n" \
  "  Certain choices of {SCALE} and {OFFSET} may result in rescaled output" \
  " values that lie outside the range {[0_1]}.  Anyway, the resulting fractional value is" \
  " clipped to {[0_1]}, then converted to an  integer sample in the range {0..OMAXVAL}, as determined by" \
  " the \"-maxval\", \"-badOut\", and \"-isMaskOut\" arguments,  and written to the output image file.\n" \
  "\n" \
  "INVALID SAMPLE VALUES\n" \
  "  The program allows a specific sample value in the input image to be declared" \
  " as meaning `undefined value'.  See" \
  " the \"-badIn\" option below.  Any output sample that depends" \
  " on `undefined' input samples will itself be `undefined'.  Pixels" \
  " outside the image domain may or" \
  " may not be considered `undefined' (see \"-replicate\" option below).  Similarly," \
  " a specific output sample value" \
  " may be used to encode `undefined' or `invalid' results.  See" \
  " the \"-badOut\" option below.\n" \
  "\n" \
  "OPERATORS\n" \
  image_window_op_INFO_OPERATORS "\n" \
  "\n" \
  "OPERATOR RANGES\n" \
  image_window_op_INFO_OP_RANGES \
  "\n" \
  "OPTIONS\n" \
  "  -op {OPKIND}.. \n" \
  "    Specifies the operator to use.  The valid alternatives are explained" \
  " in the \"OPERATORS\" section above." \
  "\n" \
  "  -smmothed\n" \
  "    This option specifies that the \"smoothed\" version of the operator" \
  " should be used.  For instance, the smoothed derivatives \"fx\" and \"fy\" are" \
  " the components of the Sobel gradient operator.  This option has no effect on" \
  " the operators \"fxy\", \"orthicity\", \"elongation\", \"average\", and \"deviation\".\n" \
  "\n" \
  "  -squared\n" \
  "    This option replaces the output of the operator by its square.\n" \
  "\n" \
  "  -scale {SCALE}\n" \
  "    This option specifies that the result of the" \
  " operator (possibly squared) be multiplied by {SCALE} before" \
  " mapping from the natural range to {[0_1]}.  The default is \"-scale 1\", that" \
  " is, no extra scaling.\n" \
  "\n" \
  "  -offset {OFFSET}\n" \
  "    This options specifies an {OFFSET} to be added to the" \
  " operator's result, after multiplying by {SCALE} and mapping the" \
  " natural range to {[0_1]}.  It may be negative.  The default is zero.\n" \
  "\n" \
  "  -maxval {OMAXVAL} \n" \
  "    Defines the nominal maximum sample value for the output" \
  " image.  Must be an integer between 1" \
  " and " stringify(PNM_MAX_MAXVAL) ".  If not specified," \
  " defaults to the input maximum sample value {IMAXVAL}.\n" \
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
  "  -keepBad \n" \
  "    This option forces an invalid output sample value" \
  " whenever the corresponding input sample is invalid.  It is effective" \
  " only if \"-badIn\" is specified, and only for some operators, such" \
  " as \"fx\", \"fy\", \"fxy\", and \"gradient\"," \
  " whose formulas do not use the center sample in the window.\n" \
  "\n" \
  "  -badOut {OBAD} \n" \
  "    If this option is present, and {OBAD} is in" \
  " the range {0..OMAXVAL}, the output sample value {OBAD}" \
  " is reserved to represent `undefined' or `invalid' operator results.  In" \
  " that case, valid operator results will be scaled to the range {0..OMAXVAL-1}" \
  " instead of {0..OMAXVAL}, and then results in the range {OBAD..OMAXVAL-1} will" \
  " be incremented by 1. Thus, for example, given \"-badOut 2\" and" \
  " {OMAXVAL = 6}, the valid operator results will be scaled so as to range" \
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
  " Valid operator results will be scaled to the" \
  " range {0..OMAXVAL}, and  invalid results will be mapped" \
  " to {OMAXVAL/2}, without any adjustment of higher values.\n" \
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
  "  -replicate\n" \
  "    This option specifies that whenever the operator needs the value of" \
  " a pixel just outside the image domain, the nearest pixel inside the" \
  " domain should be used instead.  Othwerwise, by default, those pixels" \
  " are assumed to have `undefined' value, so the formula's will be `undefined' too.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  pnmnlfilt(1), pnmconvol(1), pgmenhance(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2010-08-19 by J. Stolfi (IC-UNICAMP), based on Jef Poskanzer Netpbm image file" \
  " formats and earlier PGM filter programs.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2010-08-19 created.\n" \
  "  2020-11-15 J.Stolfi: added the \"-smoothed\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmdiffop_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

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
#include <ix.h> 
#include <uint16_image.h>
#include <float_pnm_read_stream.h>
#include <float_pnm_write_stream.h>
#include <float_image_buffer.h>
#include <image_window_op.h>
#include <sample_conv.h>
#include <affirm.h> 
#include <argparser.h> 

#define INF INFINITY

/* DATA TYPES */

typedef struct options_t 
  { char *iname;           /* Input filename ("-" for stdin). */
    char *mname;           /* Mask image name ("-" for stdin,{NULL} if not given). */ 
    char *wname;           /* Window/weight file name ("-" for stdin). */ 
    char *oname;           /* Output filename ("-" for stdout). */
    bool_t verbose;        /* TRUE says to mumble while working. */
    /* Input/output encoding and undefined value handling: */
    uint16_t maxval;       /* Output maxval (0 if not specified). */
    uint32_t badIn;        /* Input sample value that means `undefined' ({>maxval} if none). */
    bool_t keepBad;        /* TRUE forces undef output when input is undef. */
    bool_t isMaskIn;       /* TRUE if input 0 and {maxval} are to be mapped to 0.0 and 1.0. */
    uint32_t badOut;       /* Output sample value that means `undefined' ({>maxval} if none). */
    bool_t isMaskOut;      /* TRUE if output 0 and {maxval} are to mean 0.0 and 1,0. */
    bool_t replicate;      /* TRUE specifies edge replication. */
    /* Operator parameters (may be {NAN} if not applicable): */
    image_window_op_t op;  /* Kind of operator to apply. */
    bool_t smoothed;       /* TRUE to use the smoothed version of the operator. */
    bool_t squared;        /* TRUE to square the operator's result. */
    double scale;          /* Scaling factor for the result, or 1 if not given. */
    double offset;         /* Offset for the result, or 0 if not given. */
  } options_t;
  /* Arguments parsed from the command line. */

#define MAX_CHNS 3
  /* Max channels in a PNM file. */

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char* argv[]);

options_t *parse_options(int32_t argc, char **argv);

void filter_image_file
  ( FILE *ifile,
    FILE *ofile,
    options_t *o
  );
  /* Reads a PBM, PGM, or PPM image from file {ifile}, applies the operator described
    by {o}, and writes the result to {ofile}. */

void grab_input_samples
  ( int32_t c,                      /* Channel of input image. */
    int32_t x,                      /* Column of center pixel in input image. */
    int32_t y,                      /* Row of center pixel in input image. */
    float_image_buffer_t *ibuf, /* Circular row buffer of input image. */
    bool_t replicate,           /* If TRUE, extends image by edge replication. */
    int32_t wrx,                    /* Half-width of window. */
    int32_t wry,                    /* Half-height of window. */
    double wsmp[]                /* OUT: {nwx*nwy} array of pixel values in window. */
  );
  /* Stores into {wsmp[0..nw-1]} the samples taken from channel {c} of
    the input image, within the window centered at {x,y}. The window
    dimensions are assumed to be {2*wrx+1} columns by {2*wry+1}
    rows, and it is linearized by rows. If {replicate} is FALSE,
    window samples that fall outside the image's domain are set to
    {NAN}; otherwise they are copied from the nearest edge samples.
    
    Assumes that the {nwy} rows of the input image that are spanned by the
    window are stored in {buf}. */

void update_sample_range(double val, double *vminP, double *vmaxP, int32_t *nbadP);
  /* Expands the range {[*vminP _ *vmaxP]} as needed to include the
    value {val}.   If {val} is {NAN}, incrementa {*nbadP} instead. */

/* ROUTINES */

int32_t main(int32_t argc, char* argv[])
  { /* Command line arguments: */
    options_t *o = parse_options(argc, argv);
    
    /* Consistency check: */
    assert((int32_t)PNM_MAX_SAMPLE < INT_MAX - 1);

    /* Read input header and get image attributes: */
    FILE *ifile = open_read(o->iname, o->verbose);
    
    /* Write output file header: */
    FILE *ofile = open_write(o->oname, o->verbose);

    filter_image_file(ifile, ofile, o);
    
    if (ifile != stdin) { fclose(ifile); }
    if ((ofile != stdout) && (ofile != stderr)) { fclose(ofile); }
    return 0;
  }
  
void filter_image_file
  ( FILE *ifile,
    FILE *ofile,
    options_t *o
  )
  { 
    /* Choose the window size: */
    int32_t nwx = 0; /* Window width */
    int32_t nwy = 0; /* Window height */
    int32_t iwctr = 0; /* Index of window's central sample. */
    image_window_op_get_window_size(o->op, o->smoothed, &nwx, &nwy, &iwctr);
    assert(nwx % 2 == 1);
    assert(nwy % 2 == 1);
    int32_t wrx = nwx/2, wry = nwy/2; /* Half-dimensions of window. */
    int32_t wnsmp = nwx*nwy; /* Image samples per window per channel. */
    
    if (o->verbose)
      { fprintf(stderr, "using a window with %d columns and %d rows\n", nwx, nwy); }
    
    /* Read the input file header, create a buffered read stream {istr} for it: */
    float_pnm_stream_t *istr = float_pnm_read_stream_new(ifile, o->isMaskIn, o->badIn, nwy);
    int32_t rows = istr->rows;
    int32_t cols = istr->cols;
    int32_t chns = istr->chns;
    uint16_t imaxval = istr->maxval;  /* Input image's {maxval}. */
    assert(imaxval > 0);
    assert((chns > 0) && (chns <= MAX_CHNS));
    float_image_buffer_t *ibuf = istr->buf;

    if (o->verbose)
      { fprintf(stderr, "image has %d channels %d columns %d rows\n", chns, cols, rows);
        fprintf(stderr, "input maxval = %d\n", imaxval);
        if (o->badIn <= imaxval) 
          { fprintf(stderr, "input samples with value %d assumed invalid\n", o->badIn); }
      }
    
    /* Get the natural range {[vlo_vhi]} of the oeprator: */
    double vlo, vhi;
    image_window_op_get_range(o->op, o->smoothed, o->squared, &vlo, &vhi);
    if (o->verbose)
      { fprintf(stderr, "nominal operator result range = [%.6f _ %.6f]\n", vlo, vhi); }
    
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
    
    if (o->verbose)
      { fprintf(stderr, "output maxval = %d\n", omaxval);
        if (o->badOut <= imaxval) 
          { fprintf(stderr, "output samples with value %d represent invalid values\n", o->badOut); }
      }

    /* Allocate the output buffer for a single-row of quantized samples: */
    uint16_t *osmp = uint16_image_alloc_pixel_row(cols, chns);
    
    /* Actual sample ranges and {NAN} counts: */
    double vmin_in = +INF, vmax_in = -INF; int32_t nbad_in = 0; /* Among all window samples. */
    double vmin_op = +INF, vmax_op = -INF; int32_t nbad_op = 0; /* Computed op values, before scaling etc. */
    double vmin_ot = +INF, vmax_ot = -INF; int32_t nbad_ot = 0; /* Scaled/shifted op values. */
    int32_t nbad_in_ctr = 0; /* {NAN} count on central window sample only. */
    
    /* Loop on output image rows: */
    double wsmp[wnsmp];  /* Linearized window sample array, single-channel. */
    for (int32_t y = 0; y < rows; y++)
      { if (o->verbose) { fputc('!', stderr); }
        int32_t yimin = y - wry; /* Min input {y} needed to compute this row. */
        if (yimin < 0) { yimin = 0; }
        int32_t yimax = y + wry; /* Max input {y} needed to compute this row. */
        if (yimax >= rows) { yimax = rows-1; }
        
        /* Make sure that rows {yimin..yimax} are in the input buffer: */
        (void)float_pnm_read_stream_get_row(ifile, istr, yimax);
        assert(ibuf->yini <= yimin);
        assert(ibuf->ylim > yimax);
        
        /* Compute row {y} of output image: */
        for (int32_t x = 0; x < cols; x++)
          { /* Now loop on channels: */
            for (int32_t c = 0; c < chns; c++)
              { /* Get input pixel values and mask (if any) for the new window placement: */
                grab_input_samples(c, x, y, ibuf, o->replicate, wrx, wry, wsmp);
                
                if (o->verbose)
                  { for (int32_t ix = 0; ix < wnsmp; ix++)
                      { update_sample_range(wsmp[ix], &vmin_in, &vmax_in, &nbad_in); }
                  }
              
                if (o->verbose && isnan(wsmp[iwctr])) { nbad_in_ctr++; }
                
                /* Compute the floated output pixel value: */
                double oval;
                if ((o->keepBad) && (isnan(wsmp[iwctr])))
                  { oval = NAN; }
                else
                  { /* Compute the operator's value: */
                    oval = image_window_op_apply
                      ( o->op, o->smoothed, o->squared, iwctr, nwx, wsmp );
                    if (! isnan(oval))
                      { assert(isfinite(oval));
                        /* Update the actual range, for debug purposes: */
                        if (o->verbose) { update_sample_range(oval, &vmin_op, &vmax_op, &nbad_op); }

                        /* Apply the user scale: */
                        if (o->scale != 1) { oval = oval*o->scale; }

                        /* Map the natural range to {[0_1]}: */
                        oval = (oval - vlo)/(vhi - vlo);

                        /* Apply the user-requested offset: */
                        oval = oval + o->offset;
                      }
                  }
                  
                if (o->verbose) { update_sample_range(oval, &vmin_ot, &vmax_ot, &nbad_ot); }
                
                /* Quantize and store the output pixel value: */
                uint16_t oqts = pnm_quantize(oval, omaxval, o->isMaskOut, o->badOut);
                assert(oqts <= omaxval);
                osmp[x*chns + c] = oqts;
              }
          }
        /* Write row {y} of the output image: */
        pnm_write_pixels(ofile, osmp, cols, chns, omaxval, oraw, obits);
      } 
    fflush(ofile);
    if (o->verbose)
      { fputc('\n', stderr);
        fprintf(stderr, "input sample value range = [ %.8f _ %.8f ]\n", vmin_in, vmax_in);
        fprintf(stderr, "found %d bad samples in input windows (%d in central sample)\n", nbad_in, nbad_in_ctr);
        fprintf(stderr, "got %d bad operator values\n", nbad_op);
        fprintf(stderr, "raw operator value range = [ %.8f _ %.8f ]\n", vmin_op, vmax_op);
        fprintf(stderr, "got %d bad output values\n", nbad_ot);
        fprintf(stderr, "rescaled operator value range = [ %.8f _ %.8f ]\n", vmin_ot, vmax_ot);
      }

    /* Release working storage: */
    float_pnm_stream_free(istr);
    free(osmp);
  }

void update_sample_range(double val, double *vminP, double *vmaxP, int32_t *nbadP)
  {
    if (isnan(val))
      { (*nbadP)++; }
    else
      { (*vminP) = fmin((*vminP), val);
        (*vmaxP) = fmax((*vmaxP), val);
      }
  }

void grab_input_samples
  ( int32_t c,                  /* Channel of input image. */
    int32_t x,                  /* Column of center pixel in input image. */
    int32_t y,                  /* Row of center pixel in input image. */
    float_image_buffer_t *ibuf, /* Circular row buffer of input image. */
    bool_t replicate,           /* If TRUE, extends image by edge replication. */
    int32_t wrx,                /* Half-width of window. */
    int32_t wry,                /* Half-height of window. */
    double wsmp[]               /* OUT: {nwx*nwy} array of pixel values in window. */
  )
  { int32_t NC = ibuf->sz[0];
    int32_t NX = ibuf->sz[1];
    int32_t NY = ibuf->sz[2];
    
    /* Scan the window pixels: */
    int32_t kw = 0;
    for (int32_t yw = -wry; yw <= wry; yw++)
      { for (int32_t xw = -wrx; xw <= wrx; xw++)
          { /* Compute the sample's indices {xp,yp} in input img, or -1 if non-existant: */
            int32_t xp = x + xw;
            int32_t yp = y + yw;
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
                wsmp[kw] = NAN;
              }
            else
              { /* Pixel exists, get sample values from input image: */
                double *ismp = float_image_buffer_get_row(ibuf, yp);
                wsmp[kw] = ismp[xp*NC + c];
              }
            kw++;
          }
      }
  }

/* AUXILIARY PROCS */

options_t *parse_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
    
    /* Get keyword arguments: */

    /* Operator kind and parameters: */
    if (argparser_keyword_present(pp, "-op"))
      { 
        char *chop = argparser_get_next_non_keyword(pp); 
        o->op = image_window_op_from_string(chop);
        if (o->op >= image_window_op_NUM)
          { argparser_error(pp, "invalid operator kind"); }
      }
    else
      { argparser_error(pp, "missing \"-op\" option"); }

   o->smoothed = argparser_keyword_present(pp, "-smoothed");

   o->squared = argparser_keyword_present(pp, "-squared");
    
    if (argparser_keyword_present(pp, "-scale"))
      { o->scale = argparser_get_next_double(pp, -1.0e+20, +1.0e+20); }
    else
      { o->scale = 1.0; }
    
    if (argparser_keyword_present(pp, "-offset"))
      { o->offset = argparser_get_next_double(pp, -1.0e+20, +1.0e+20); }
    else
      { o->offset = 0.0; }
    
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (uint16_t)argparser_get_next_int(pp, 0, PNM_FILE_MAX_MAXVAL); }
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
    
    o->verbose = argparser_keyword_present(pp, "-verbose");
    
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
