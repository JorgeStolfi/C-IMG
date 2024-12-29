#define PROG_NAME "pnmadjust"
#define PROG_DESC "color correction for digital and scanned photographs"
#define PROG_VERS "1.0"

#define pnmadjust_C_COPYRIGHT "Copyright © 2003 by the State University of Campinas (UNICAMP)"

/* Last edited on 2024-12-21 12:00:05 by stolfi */

/* TO DO: !!! Unify the documentation with that of {pnmfield.c} !!! */

#define PROG_HELP \
  PROG_NAME "\\\n" \
  "    [ -inGamma {TRIPLET} ].. \\\n" \
  "    [ -outGamma {TRIPLET} ].. \\\n" \
  "    [ -kappa {TRIPLET} ].. \\\n" \
  "    [ -black {COLOR} | -varBlack {FIELDSPECS} ] \\\n" \
  "    [ -white {COLOR} | -varWhite {FIELDSPECS} ] \\\n" \
  "    [ -showWhite | -showBlack ] \\\n" \
  "    [ -debug {DEBX} {DEBY} ] \\\n" \
  "    [ {PNMFILE} ] \\\n" \
  "  Each {TRIPLET} is three real numbers.\n" \
  "  Each {COLOR} is\n" \
  "    " frgb_parse_color_HELP "\n" \
  "  See the pnmfield(1) documentation for the {FIELDSPECS} syntax.\n"

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
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "EXAMPLES\n" \
  "\n" \
  PROG_INFO_EXMP "\n" \
  "SEE ALSO\n" \
  "  pnmgamma(1), pgmnorm(1), pnm(5)." \
  "AUTHOR\n" \
  "  Created in 2002 by Jorge Stolfi, IC-UNICAMP." \
  "  Based on ColorCorrect.m3 by Jorge Stolfi, jun/1995." \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2006-11-21 by J. Stolfi, IC-UNICAMP.  Folded the manpage into the program.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmadjust_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "   Reads an input PGM/PPM file containing a digital" \
  " or scanned photograph.  Applies some commonly needed color" \
  " corrections to each pixel, and writes the result to standard" \
  " output.\n" \
  "\n" \
  "  If the input is color (PPM) or grayscale (PGM), the output" \
  " will be of the same type and depth.\n" \
  "\n" \
  "  The program assumes that the input is an `ideal' image of" \
  " a real scene that has been modified by the following processes," \
  " in order: (1) multiplication by a `light field' image; (2) addition" \
  " of a `dark field' image; (3) non-linear (`gamma') encoding; (4)" \
  " quantization.  The program undoes these perturbations, in the" \
  " reverse order, to obtain the ideal intensities in a linear energy" \
  " scale. Then it applies optional brightness (`kappa') and saturation" \
  " enhancements. Finally, it applies another `gamma' encoding and" \
  " re-quantizes the result.  The gamma decoding and encoding is" \
  " performed by {sample_conv_gamma} (q.v.) with" \
  " user-given {expo} parameters and a standard {bias}."
  
#define PROG_INFO_OPTS \
  "  In what follows, the {TRIPLET} parameters are three real numbers," \
  " one for each color channel (red, green, blue).  The {COLOR} parameters" \
  " are " frgb_parse_color_INFO "\n" \
  "\n" \
  "  The {COLOR} components refer to sample values in the" \
  " input image, only scaled to the [0 _ 1] range; they are thus assumed to" \
  " be affected by light-field factors, dark-level noise, and gamma" \
  " encoding.\n" \
  "\n" \
  "  -inGamma {TRIPLET}\n" \
  "    Specifies the exponent of the nonlinear encoding (`gamma')" \
  " for each channel of the input image.\n" \
  "\n" \
  "  -black {COLOR}\n" \
  "    Specifies the pixel value in the actual image which should be" \
  " black (RGB = 0,0,0) in the ideal image.  Use this parameter to" \
  " compensate for uniform `dark field' bias, atmospheric glare, etc..\n" \
  "\n" \
  "  -varBlack {FIELDSPECS}\n" \
  "    Specifies that the `black' reference value changes from place" \
  " to place in the image.  See the pnmfield(1) manpage for the{FIELDSPECS}" \
  " syntax.  Use this option to compensate for non-uniform dark-field bias," \
  " such as glare from a dirty window.  This option and \"-black\" are" \
  " mutually exclusive.\n" \
  "\n" \
  "  -white {COLOR}\n" \
  "    Specifies the pixel value in the actual image which should be" \
  " white (RGB = 1,1,1) in the ideal image. Use this parameter to compensate" \
  " for wrong exposure or a non-white light source.\n" \
  "\n" \
  "  -varWhite {FIELDSPECS}\n" \
  "    Specifies that the `white' reference value changes from place to" \
  " place in the image. See the the pnmfield(1) manpage for the {FIELDSPECS}" \
  " syntax.  Use this option to compensate for non-uniform lighting of the" \
  " scene; in particular, to remove the yellow/blue stripes caused by" \
  " fluorescent lighting.  This option and \"-white\" are mutually exclusive.\n" \
  "\n" \
  "  -kappa {TRIPLET}\n" \
  "    Specifies a nonlinear brightness adjustment (`kappa correction')" \
  " to be applied to the ideal image, after compensating for input" \
  " gamma-encoding, black-field level, and white-field level" \
  " (uniform or variable). In general, the minimum (0) and maximum (1)" \
  " intensity values remain unchanged, while the middle ones are" \
  " increased ({kappa > 1}) or reduced ({kappa < 1}).  If the same" \
  " {kappa} is specified for all three channels, the original hue" \
  " and saturation are preserved.\n" \
  "\n" \
  "  -saturation {NUMBER}\n" \
  "    Specifies that the saturation of colors the ideal image" \
  " should be increased ({NUMBER > 1}) or reduced ({NUMBER < 1})" \
  " without changing their intensities.\n" \
  "\n" \
  "  -outGamma {TRIPLET}\n" \
  "    Specifies the nonlinear exponent to be used when encoding" \
  " each channel of the the corrected and adjusted intensities before writing them" \
  " to the output file. See the \"-inGamma\" option for details.\n" \
  "\n" \
  "  -showWhite\n" \
  "  -showBlack\n" \
  "    These mutually exclusive options specify that the output should" \
  " be a diagnostic image showing the white reference field or the black" \
  " reference field, instead of the adjusted image. Note that the field" \
  " is clipped to the [0_1] range for image output.\n" \
  "\n" \
  "  -debug {DEBX} {DEBY}\n" \
  "    If present, this option requests debugging information for the" \
  " pixel in column {DEBX} from left and row {DEBY} from top (both" \
  " counted from 0)."
  
#define PROG_INFO_EXMP \
  "  The following example corrects a linear-scale image" \
  " \"foo.ppm\" for a uniform dark-field bias and a light" \
  " field that varies according to two orthogonal waves," \
  " with periods (300,100) and (-50,150).  The resulting" \
  " image is re-encoded with gamma exponent = 2.2 and written" \
  " out to file\"bar.ppm\":\n" \
  "\n" \
  "    pnmadjust \\\n" \
  "      -inGamma      1.0 1.0 1.0       \\\n" \
  "      -black        030 025 020 / 255 \\\n" \
  "      -varWhite wavePair  \\\n" \
  "         0120 0320  240 235 230 / 255 \\\n" \
  "         0420 0420  201 198 180 / 255 \\\n" \
  "         0070 0470  185 180 178 / 255 \\\n" \
  "      -outGamma     2.2 2.2 2.2       \\\n" \
  "      < foo.ppm \\\n" \
  "      > bar.ppm\n" \
  "\n" \
  "  The following example increases the brightness by 50%" \
  " and the saturation by 10%. The input and output images" \
  " have gamma = 2.2:\n" \
  " \n" \
  "    pnmadjust \\\n" \
  "      -inGamma    2.2 2.2 2.2 \\\n" \
  "      -kappa      1.5 1.5 1.5 \\\n" \
  "      -saturation 1.1         \\\n" \
  "      -outGamma   2.2 2.2 2.2 \\\n" \
  "      < foo.ppm \\\n" \
  "      > bar.ppm\n"

#include <stdio.h>
#include <values.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <argparser.h>
#include <bool.h>
#include <frgb.h>
#include <jsfile.h>
#include <frgb_ops.h>
#include <colorfield.h>
#include <sample_conv.h>
#include <jspnm.h>
#include <uint16_image.h>

/* Maximum rows, columns, elements (to avoid absurd allocs): */
#define MAX_PIXELS (16*256*256*256)   
#define MAX_SIZE MAX_PIXELS

#define BT_BIAS (sample_conv_BT709_BIAS) 
/* A {bias} parameter that approximates BT.709 when
  used with expo {0.45} or {1/0.45}. */

/* COMMAND LINE ARGUMENTS */

typedef struct options_t 
  { char *f_in_name;       /* Input filename, or "-" for stdin. */
    frgb_t inGamma_expo;        /* Gamma of input image, per channel. */
    cfld_args_t *black;    /* The {black} color field, or NULL. */
    cfld_args_t *white;    /* The {white} color field, or NULL. */
    double globalKappa;    /* Kappa correction for intensity. */
    frgb_t channelKappa;   /* Additional kappa correction, per channel. */
    double saturation;     /* Saturation enhancement factor. */
    frgb_t outGamma_expo;       /* Gamma of output image, per channel. */
    bool_t showWhite;      /* TRUE to show the white reference field. */
    bool_t showBlack;      /* TRUE to show the black reference field. */
    int debug_col;         /* Column of pixel to debug, or -1 if none. */
    int debug_row;         /* Row of pixel to debug, or -1 if none. */
  } options_t;

/* INTERNAL PROTOTYPES */

bool_t eqn(float a[], float b[], int n);
  /* TRUE iff {a[i]=b[i]} for {i=0..n-1}. */

int main(int argc, char **argv);

options_t *parse_options(int argc, char **argv);

int scale_is_trivial(options_t *o, int chns);
  /* TRUE if the input requires no black- or white-level correction
    other than input gamma-correction. */
    
cfld_params_t *compute_color_field
  ( cfld_args_t *cfargs, 
    frgb_t *inGamma_expo, 
    int chns
  );
  /* Computes the parameters of a gamma-corrected reference color
    field. The result depends on {cfargs} and {inGamma_expo}.
    If {chns == 1}, all color arguments are reduced to grayscale. */

cfld_params_t *compute_white_field
  ( cfld_args_t *wfargs, 
    frgb_t *inGamma_expo, 
    cfld_params_t *black, 
    int chns
  );
  /* Computes the parameters of the gamma- and black-corrected `white'
    reference field. The result depends on the white firls arguments 
    {wfargs}, {inGamma_expo},  and the previously computed black-field map {black}.
    If {chns == 1}, all color arguments are reduced to grayscale. */

void read_input_file_header(FILE *rd, uint16_image_t **hd, pnm_format_t *fmt);
  /* Reads the header of a PBM/PGM/PPM file, returns the data in {hd}
    and the file format in {*fmt}. */
  
void make_output_file_header
  ( uint16_image_t *hd_in, 
    pnm_format_t fmt_in, 
    uint16_image_t **hd_ot, 
    pnm_format_t *fmt_ot
  );
  /* Creates the header {hd_ot} of the output file, given that {hd_in} 
    of the input file. */

void write_output_file_header(FILE *wr, uint16_image_t *hd_ot, pnm_format_t fmt_ot);
  /* Writes the header {hd_ot} to the output file {wr}, with format {fmt_ot}
    and {forceplain = FALSE}. */
  
void process_pixels(FILE *rd, options_t *o, FILE *wr);
  /* Reads body of input image from file {rd}, computes output pixels,
    writes them to {wr}. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  { 
    /* Parse command line arguments: */
    options_t *o = parse_options(argc, argv);

    /* Open files: */
    FILE *rd = open_read(o->f_in_name, TRUE);
    fprintf(stderr, "writing result to (stdout)...\n");
    FILE *wr = stdout;
    
    /* Process image: */
    process_pixels(rd, o, wr);

    if (wr == stdout) { fflush(wr); } else { fclose(wr); }
    if (rd != stdin) { fclose(rd); }
    fprintf(stderr, "done.\n");
    return 0;
  }

void process_pixels(FILE *rd, options_t *o, FILE *wr)
  { 
    /* Read the input file header: */
    int rows, cols, chns;
    uint16_t imaxval;
    bool_t iraw, ibits;
    pnm_format_t iformat;
    pnm_read_header(rd, &cols, &rows, &chns, &imaxval, &iraw, &ibits, &iformat);
    assert((chns == 1) || (chns == 3));
    
    /* Choose the output maxval and format, and write the output file header: */
    uint16_t omaxval = imaxval;
    pnm_format_t oformat;
    bool_t oraw, obits;
    pnm_choose_output_format(omaxval, chns, (! iraw), &oformat, &oraw, &obits);
    pnm_write_header(wr, cols, rows, omaxval, oformat);
    
    /* Input and output scaling parameters: */
    double zero_in = 0;
    double unit_in = (double)imaxval;
    double scale_in = unit_in - zero_in;
    
    double zero_ot = 0;
    double unit_ot = (double)omaxval;
    double scale_ot = unit_ot - zero_ot;
    
    /* Required pixel processing stages: */
    bool_t do_in_gamma = ! eqn(o->inGamma_expo.c, frgb_Ones.c, chns);
    bool_t do_scale = (! scale_is_trivial(o, chns)) | o->showWhite | o->showBlack;
    bool_t do_chan_kappa = ! eqn(o->channelKappa.c, frgb_Ones.c, chns);
    bool_t do_glob_kappa = ! (o->globalKappa == 1.0);
    bool_t do_saturation = ! (o->saturation == 1.0);
    bool_t do_ot_gamma = ! eqn(o->outGamma_expo.c, frgb_Ones.c, chns);

    frgb_t *inGamma_expo = (do_in_gamma ? &(o->inGamma_expo) : NULL);

    /* Black and white reference fields: */
    cfld_params_t *blackField = compute_color_field(o->black, inGamma_expo, chns);
    cfld_params_t *whiteField = compute_color_field(o->white, inGamma_expo, chns);
    
    /* Alocate row buffers: */
    uint16_t *smp_row_in = uint16_image_alloc_pixel_row(cols, chns);
    uint16_t *smp_row_ot = uint16_image_alloc_pixel_row(cols, chns);

    /* Loop on rows: */
    uint16_t *sp_in, *sp_ot;
    int row, col, i;
    frgb_t fv;
    int iv[chns];
    for (row = 0; row < rows; ++row)
      { pnm_read_pixels(rd, smp_row_in, cols, chns, imaxval, iraw, ibits);
        sp_in = &(smp_row_in[0]); sp_ot = &(smp_row_ot[0]);
        for (col = 0; col < cols; ++col)
          { /* Shall we debug this pixel? */
            bool_t debug = (col == o->debug_col) & (row == o->debug_row);
            frgb_DEBUG = debug;
            for (i = 0; i < chns; i++)
              { fv.c[i] = frgb_floatize(*sp_in, imaxval, zero_in, scale_in); sp_in++; }
            frgb_debug("ip", col, row, &fv, chns, "\n");
            /* Get the local black and white values, if needed: */
            frgb_t locBlack, locWhite; /* Local black and white values. */
            if (do_scale)
              { /* Compute the local black and white reference: */
                if (blackField != NULL) 
                  { cfld_eval(blackField, col, row, &locBlack, chns); }
                else
                  { locBlack = frgb_Black; }
                if (whiteField != NULL) 
                  { cfld_eval(whiteField, col, row, &locWhite, chns); }
                else
                  { locWhite = frgb_White; }
                frgb_debug("lb", col, row, &locBlack, chns, "\n");
                frgb_debug("lw", col, row, &locWhite, chns, "\n");
              }
            for (i = 0; i < chns; i++)
              { double fvc = fv.c[i];
                if (do_in_gamma)
                  { /* Convert input values to linear scale: */
                    fvc = sample_conv_gamma(fvc, o->inGamma_expo.c[i], BT_BIAS);
                  }
                if (do_scale)
                  { /* Apply linear black- and white-level correction: */
                    double blc = locBlack.c[i];
                    double wht = locWhite.c[i];
                    double dif = wht - blc;
                    if (o->showWhite)
                      { fvc = wht; }
                    else if (o->showBlack)
                      { fvc = blc; }
                    else if (dif <= 1.0e-6)
                      { fvc = 0.5; }
                    else 
                      { fvc = (fvc <= blc ? 0.0 : fvc - blc); 
                        fvc = (fvc >= dif ? 1.0 : fvc / dif);
                      }
                  }
                if (do_chan_kappa)
                  {/* Apply `kappa' adjustment to each channel: */
                    fvc = frgb_apply_kappa_gray(fvc, o->channelKappa.c[i]);
                  }
                fv.c[i] = fvc;
              }
            frgb_debug("cp", col, row, &fv, chns, "\n");
            if (do_scale || do_glob_kappa || do_saturation)
              { /* Apply global `kappa' and saturation adjustments, and clip: */
                frgb_apply_glob_kappa_sat_clip(&fv, o->globalKappa, o->saturation);
              }
            frgb_debug("ap", col, row, &fv, chns, "\n");
            if (do_ot_gamma)
              { for (i = 0; i < chns; i++)
                  { /* Convert linear scale to output pixels: */
                    fv.c[i] = sample_conv_gamma(fv.c[i], 1/o->outGamma_expo.c[i], BT_BIAS);
                  }
              }
            frgb_debug("op", col, row, &fv, chns, "\n");
            for (i = 0; i < chns; i++) 
              { iv[i] = frgb_quantize(fv.c[i], zero_ot, scale_ot, omaxval); }
            frgb_debug_int_pixel("oq", col, row, iv, chns, "\n");
            for (i = 0; i < chns; i++) 
              { (*sp_ot) = iv[i]; sp_ot++; }
            if (frgb_DEBUG) { fprintf(stderr, "\n\n"); }
          }
        pnm_write_pixels(wr, smp_row_ot, cols, chns, omaxval, oraw, obits);
      }
  }
  
int scale_is_trivial(options_t *o, int chns)
  { return (o->black == NULL) && (o->white == NULL); }

cfld_params_t *compute_color_field
  ( cfld_args_t *cfargs, 
    frgb_t *inGamma_expo, 
    int chns
  )
  { if (cfargs == NULL)
      { return NULL; }
    else
      { auto frgb_t adjust_color_arg(frgb_t *v, int col, int row);

        frgb_t adjust_color_arg(frgb_t *v, int col, int row)
          { return frgb_correct_arg(v, inGamma_expo, (chns == 1)); }

        cfld_params_t *bfp = cfld_compute_params(cfargs, adjust_color_arg, FALSE);
        return bfp;
      }
  }

cfld_params_t *compute_white_field
  ( cfld_args_t *wfargs, 
    frgb_t *inGamma_expo, 
    cfld_params_t *black, 
    int chns
  )
  { if (wfargs == NULL)
      { return NULL; }
    else
      { auto frgb_t adjust_white_arg(frgb_t *v, int col, int row);

        frgb_t adjust_white_arg(frgb_t *v, int col, int row)
          { frgb_t fv = frgb_correct_arg(v, inGamma_expo, (chns == 1));
            if (black != NULL)
              { frgb_t bv;
                cfld_eval(black, col, row,  &bv, chns);
                int i;
                for (i = 0; i < 3; i++) { fv.c[i] -= bv.c[i]; }
              }
            return fv;
          }

        cfld_params_t *wfp = cfld_compute_params(wfargs, adjust_white_arg, TRUE);
        return wfp;
      }
  }

options_t *parse_options(int argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 

    /* Parse command line options: */
    if (argparser_keyword_present(pp, "-inGamma"))
      { o->inGamma_expo = frgb_parse(pp, 0.1, 10.0); }
    else
      { o->inGamma_expo = (frgb_t){{1.0, 1.0, 1.0}}; }
    
    if (argparser_keyword_present(pp, "-black"))
      { frgb_t color = frgb_parse_color(pp);
        if ((color.c[0] == 0) && (color.c[1] == 0) && (color.c[2] == 0))
          { o->black = NULL; }
        else
          { o->black = cfld_make_args_uniform(&color); }
      }
    else if (argparser_keyword_present(pp, "-varBlack"))
      { o->black = cfld_parse(pp);  }
    
    if (argparser_keyword_present(pp, "-white"))
      { frgb_t color = frgb_parse_color(pp);
        if ((color.c[0] == 1) && (color.c[1] == 1) && (color.c[2] == 1))
          { o->white = NULL; }
        else
          { o->white = cfld_make_args_uniform(&color); }
      }
    else if (argparser_keyword_present(pp, "-varWhite"))
      { o->white = cfld_parse(pp); }
    
    if (argparser_keyword_present(pp, "-kappa"))
      { o->channelKappa = frgb_parse(pp, 0.01, 100.0);
        frgb_t *ck = &(o->channelKappa);
        double avg = pow(ck->c[0]*ck->c[1]*ck->c[2], 1.0/3.0);
        o->globalKappa = avg;
        int i;
        for (i = 0; i < 3; i++) { ck->c[i] = ck->c[i]/avg; }
      }
    else
      { o->channelKappa = (frgb_t){{1.0, 1.0, 1.0}};
        o->globalKappa = 1.0;
      }
    
    if (argparser_keyword_present(pp, "-saturation"))
      { o->saturation = argparser_get_next_double(pp, 0.0, 100.0); }
    else
      { o->saturation = 1.0; }
    
    if (argparser_keyword_present(pp, "-outGamma"))
      { o->outGamma_expo = frgb_parse(pp, 0.1, 10.0); } 
    else
      { o->outGamma_expo = (frgb_t){{1.0, 1.0, 1.0}}; }

    if (argparser_keyword_present(pp, "-showWhite"))
      { o->showWhite = TRUE; o->showBlack = FALSE; } 
    else if (argparser_keyword_present(pp, "-showBlack"))
      { o->showBlack = TRUE; o->showWhite = FALSE; } 
    else
      { o->showBlack = o->showWhite = FALSE; }

    if (argparser_keyword_present(pp, "-debug"))
      { o->debug_col = argparser_get_next_int(pp, -1, MAX_SIZE+1);
        o->debug_row = argparser_get_next_int(pp, -1, MAX_SIZE+1);
      } 
    else
      { o->debug_col = o->debug_row = -1 ; }

    argparser_skip_parsed(pp);

    if (argparser_next(pp) != NULL) 
      { o->f_in_name = argparser_get_next(pp); }
    else
      { o->f_in_name = "-"; }

    argparser_finish(pp);
        
    return o;
  }

bool_t eqn(float a[], float b[], int n)
  { int i; 
    for (i = 0; i < n; i++) { if (a[i] != b[i]) { return FALSE; } }
    return TRUE;
  }
