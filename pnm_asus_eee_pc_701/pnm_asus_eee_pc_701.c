#define PROG_NAME "pnm_asus_eee_pc_701"
#define PROG_DESC "Reformat Asus Eee-Pc 701 camera images "
#define PROG_VERS "1.0"

/* !!! Finish !!! */

#define pnm_asus_eee_pc_701_C_COPYRIGHT "Copyright © 2008 by the State University of Campinas (UNICAMP)"
/* Last edited on 2024-12-21 12:00:20 by stolfi */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -babla \\\n" \
  "  [ -verbose ]"

#define PROG_INFO_DESC \
  "  The program reads a raw PNM image captured with the Asus Eee-PC 701 camera.  It outputs a nicer one.\n" \
  "\n" \
  "OUTPUT\n" \
  "  Bla bla bla..."

#define PROG_INFO_OPTS \
  "  -verbose\n" \
  "    If this optional keyword is present, the program prints" \
  " various diagnostic messages."
  
#define PROG_INFO_EXMP \
  "  The Asus Eee-Pc 701 camera has two program-settable" \
  " parametes,  brightness {bbb} (an integer from 0 to 255) and" \
  " contrast {ccc} (an integer from 1 to 31).   The format of" \
  " the captured frames is YUV 4:2:2 where each sample is an" \
  " integer in 0 to 255.  At some point the camera computes" \
  " an internal luminance {Z} in the range {16..255}, where" \
  " {Z==16} is actually some positive luminance.  Then" \
  " the output luminance {Y} is computed from {Z} by the formula\n" \
  "\n" \
  "   {Y = max(16, min(255, (Z ± G(2))*ccc/16 + (bbb - 128)))}.\n" \
  "\n" \
  " Matching this equation to the camera model above we get\n" \
  "\n" \
  "   {WHITE = 255}\n" \
  "   {BLACK = 016}\n" \
  "   {BRGHT = (bbb - 128 - BLACK)/(WHITE-BLACK)}\n" \
  "   {CTRST = ccc/16/(WHITE-BLACK)}\n" \
  "   {SIGMA ~ 2}\n" \
  "   {GAMMA = 1.0000}.\n" \
  "\n" \
  " In particular, we have\n" \
  "\n" \
  "   {bbb} {BRGHT} \n" \
  "   ----- ------- \n" \
  "    025  -0.4979 \n" \
  "    050  -0.3933 \n" \
  "    075  -0.2887 \n" \
  "    100  -0.1841 \n" \
  "    125  -0.0795 \n" \
  "    150  +0.0251 \n" \
  "    175  +0.1297 \n" \
  "    200  +0.2343 \n" \
  "    225  +0.3389 \n" \
  "    250  +0.4435 \n" \
  "\n" \
  "   {ccc} {CTRST} \n" \
  "   ----- ------- \n" \
  "    006  +0.0016 \n" \
  "    012  +0.0031 \n" \
  "    018  +0.0047 \n" \
  "    024  +0.0063 \n" \
  "    030  +0.0078 \n"

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
  "EXAMPLE: THE ASUS EEE-PC 701\n" \
  PROG_INFO_EXMP "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnmscale(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created nov/2008 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  nov/2008 Created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnm_asus_eee_pc_701_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <jsfile.h>
#include <jspnm.h>
#include <vec.h>
#include <uint16_image.h>
#include <sample_conv_hdyn.h>
#include <float_image.h>
#include <float_image_hdyn.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <argparser.h>

typedef struct options_t 
  { bool_t blabla;         /* Input file names. */
    /* Debugging options: */
    bool_t verbose;      /* TRUE to print global statistics. */
  } options_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char **argv);

options_t *get_options(int argc, char **argv);
  /* Gets the arguments from the command line and returns them as an
    {options_t} record {o}. */

uint16_image_t *read_input_image(FILE *rd, bool_t verbose);
  /* Reads a PBM/PGM image file from {rd}.
    If {verbose} is true, prints image statistics to {stderr}. */

void write_output_image
  ( FILE *wr, 
    float_image_t *fim, 
    uint16_t maxval,
    bool_t verbose
  );
  /* Writes the float image {fim} to {wr} as a PBM/PGM/PPM image file.
    Samples are converted from the range {[0 _ 1]} to {0 .. maxval-1}
    with {NaN} or {±INF} values being converted to {maxval}.
    
    If {verbose} is true, prints image statistics to {stderr}. */
    
uint16_image_t *float_image_to_uint16_image_safe
  ( float_image_t *fim, 
    int chns,
    double lo[], 
    double hi[], 
    int ch[],
    uint16_t maxval, 
    uint16_t bad,
    bool_t yup,
    bool_t verbose
  );
  /* Converts the image {fim} to a PNM image, like
    {float_image_to_uint16_image}, except that it safely handles infnite
    and NAN values, and uses a slighlty different mapping for other
    values.
    
    Namely, this procedure will map {lo[c]} to {1} and {hi[c]} to
    {maxval-1}, for each channel {c}. Values below {lo[c]} (including
    {-INF}) are mapped to 0; values above {hi[c]} (including {+INF})
    are mapped to 1; and {NAN} values are mapped to {bad}. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    /* Parse command line options: */
    options_t *o = get_options(argc, argv);
    int NI = o->fname.ne; /* Number of input images. */
    assert(o->params.ne == NI);

    /* Read input images, get/check dimensions: */
    uint16_image_t *pim[NI];
    int cols, rows;
    int k;
    for (k = 0; k < NI; k++)
      { FILE *rd = open_read(o->fname.e[k], o->verbose);
        pim[k] = read_input_image(rd, o->verbose);
        fclose(rd);
        demand(pim[k]->chns == 1, "input images must be monochromatic");
        /* Set/check {cols,rows}: */
        if (k == 0)
          { cols = pim[k]->cols; rows = pim[k]->rows; }
        else
          { demand(pim[k]->cols == cols, "inconsistent image width");
            demand(pim[k]->rows == rows, "inconsistent image height");
          }
      }
    
    /* Compute combined brightness image: */
    float_image_t *fim = float_image_hdyn_combine(NI, pim, o->params.e, o->verbose);
    
    /* Write the filtered image, scaled and quantized: */
    uint16_t maxval_ot = 65535;
    write_output_image(stdout, fim, maxval_ot, o->verbose);

    if (o->verbose) { fprintf(stderr, "done."); }
    return 0;
  }

uint16_image_t *read_input_image(FILE *rd, bool_t verbose)
  { uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    return pim;
  }

void write_output_image
  ( FILE *wr, 
    float_image_t *fim, 
    uint16_t maxval,
    bool_t verbose
  )
  { int chns = fim->sz[0];
    demand(chns == 1, "image must be monchromatic");
    double vLo[chns];
    double vHi[chns];
    int c;
    for (c = 0; c < chns; c++) { vLo[c] = 16; vHi[c] = 255; }
    uint16_image_t *pim = float_image_to_uint16_image_safe(fim, chns, vLo, vHi, NULL, maxval, 0, TRUE, verbose);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pim, forceplain, verbose);
    uint16_image_free(pim);
  }

uint16_image_t *float_image_to_uint16_image_safe
  ( float_image_t *fim, 
    int chns,
    double lo[], 
    double hi[], 
    int ch[],
    uint16_t maxval, 
    uint16_t bad,
    bool_t yup, 
    bool_t verbose
  )
  { /* Get image dimensions: */
    int NX = fim->sz[1];
    int NY = fim->sz[2];
    
    /* Channel counts: */
    int fchns = fim->sz[0]; /* Num channels in float image. */
    int ichns = chns;       /* Num channels in integer image. */
    
    /* Allocate PGM/PPM image: */
    uint16_image_t *iim = uint16_image_new(NX, NY, ichns);
    
    /* Set max sample value in integer image: */
    iim->maxval = maxval;
    
    /* Channel indexing variables: */
    int k; /* Channel of integer image. */
    int c; /* Channel of float image. */
    
    /* Input and output range registers: */
    float vmin[ichns], vmax[ichns];         /* Float pixel range. */
    sample_uint_t imin[ichns], imax[ichns]; /* Int pixel range. */
    int clo[ichns], chi[ichns];             /* Counts of lo-clipped and hi-clipped pixels. */
    for (k = 0; k < ichns; k++) 
      { clo[k] = chi[k] = 0;
        vmin[k] = +INF;
        vmax[k] = -INF; 
        imin[k] = maxval;
        imax[k] = 0;  
      }
    
    /* Convert pixels, store in {iim}, keep statistics: */
    int x, y;
    for(y = 0; y < NY; y++)
      { int ppmy = (yup ? NY - 1 - y : y);
        uint16_t *prow = iim->smp[ppmy];
        for(x = 0; x < NX; x++)
          { /* Convert float pixel {fpxy[c..c+2]} to integer pixel {ipxy[0..2]}, keep stats: */
            for (k = 0; k < ichns; k++)
              { double lok = (lo == NULL ? 0.0 : lo[k]);
                double hik = (hi == NULL ? 1.0 : hi[k]);
                c = (ch == NULL ? k : ch[k]);
                float v = ((c < 0) || (c >= fchns) ? 0.0 : float_image_get_sample(fim, c, x, y));
                sample_uint_t iv;
                if (isnan(v))
                  { iv = bad; }
                else if (v == -INF)
                  { iv = 0; }
                else if (v == +INF)
                  { iv = maxval; }
                else
                  { bool_t isMask = FALSE; /* Assume smooth distribution of sample values. */
                    iv = 1 + 
                      sample_conv_quantize
                        ( v, maxval-2, isMask, lok, hik, 
                          &(vmin[k]), &(vmax[k]), &(clo[k]), &(chi[k]), &(imin[k]), &(imax[k])
                        );
                  }
                (*prow) = iv;
                prow++;
              }
          }
      }
    
    if (verbose)
      { /* Print statistics: */
        long int NPIX = ((long int)NX)*((long int)NY);
        fprintf(stderr, "  %ld pixels in float image\n", NPIX);
        if (NPIX > 0)
          { for (k = 0; k < chns; k++)
              { double lok = (lo == NULL ? 0.0 : lo[k]);
                double hik = (hi == NULL ? 1.0 : hi[k]);
                c = (ch == NULL ? k : ch[k]);
                sample_conv_print_quantize_stats
                  ( c, k, vmin[k], vmax[k], lok, hik, clo[k], chi[k], maxval, imin[k], imax[k]);
              }
          }
      }
    
    return iim;
  }


#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

#define MAX_SIZE (32*1024)
  /* A limit on image size, to avoid humongous mallocs. */

options_t *get_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
     
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    o->fname = string_vec_new(10);
    o->params = params_vec_new(10);
    
    int ni = 0;
    while (argparser_keyword_present(pp, "-image"))
      { 
        string_vec_expand(&(o->fname), ni);
        params_vec_expand(&(o->params), ni);
        o->fname.e[ni] = argparser_get_next(pp);
        o->params.e[ni].brght = argparser_get_next_double(pp, -100000.0, +100000.0);
        o->params.e[ni].ctrst = argparser_get_next_double(pp, -100000.0, +100000.0);
        o->params.e[ni].sigma = argparser_get_next_double(pp, 0.00001, 100000.0);
        o->params.e[ni].gamma = argparser_get_next_double(pp, 0.1, 10.0);
        o->params.e[ni].black = argparser_get_next_int(pp, 0, PNM_MAX_SAMPLE-1);
        o->params.e[ni].white = argparser_get_next_int(pp, o->params.e[ni].black, PNM_MAX_SAMPLE);
        ni++;
     }
    string_vec_trim(&(o->fname), ni);
    params_vec_trim(&(o->params), ni);
    
    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }

vec_typeimpl(params_vec_t,params_vec,params_t);

