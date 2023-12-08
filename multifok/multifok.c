#define PROG_NAME "multifok"
#define PROG_DESC "Multi-focus stereo microscopy"
#define PROG_VERS "1.0"

// Last edited on 2023-11-25 18:28:47 by stolfi

#define multifok_C_COPYRIGHT \
    "© 2017 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -fileFormat {FILE_FORMAT} \\\n" \
  "    { -image {Z_k} {IMAGE_k} }.. \\\n" \
  "    [ -gammaDec {GAMMADEC} ] [ -bias {BIAS} ] \\\n" \
  "    [ -noise {NOISE} ] \\\n" \
  "    [ -spread {SPREAD} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    -outPrefix {OPREF} \\\n" \
  "    < {PAIRFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "   Reads a stack of images of an object (/frames/), taken with a camera with" \
  " limited depth of focus, all from the same direction and with same" \
  " lighting, aperture, and speed, but with different object-camera" \
  " distances.  Determines which parts of each image are in best" \
  " focus, and combines them into a single image.  As a side" \
  " effect, returns a /focus map/ for each image, and an" \
  " approximate /height map/ for the object.\n" \
  "\n" \
  "   The input files must be in the specified format, all with" \
  " the same size, and must be aligned and scaled so that the" \
  " focused parts of the object appear in the same scale in every image.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The program writes to \"{OPREF}-final.png\" an image that merges" \
  " the best parts of each frame into a single image that is hoped" \
  " to be sharp everywhere.  The image will have the same number of" \
  " channels as the input images.\n" \
  "\n" \
  "  The program also writes to \"{OPREF}-{KKK}-foc.png\" a single-channel PNG" \
  " image whose value at some pixel is a non-negative score that" \
  " indicates the apparent quality of focus around that pixel.  The" \
  " field {KKK} is the index of the image in the stack, as three decimal digits.\n" \
  "\n" \
  "  The program also writes to \"{OPREF}-{KKK}-mix.png\" a single-channel PNG" \
  " image whose value at some pixel is the mixing coefficient of each image in the final image.\n" \
  "\n" \
  "  The program also writes to \"{OPREF}-height.fni\" a three-channel image" \
  " with an approximate height map of the object.  The first channel" \
  " is the estimated mean Z-coordinate of the surface visible in" \
  " each pixel.  The second channel is the estimated absolute" \
  " deviation of the height.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -fileFormat {FILE_FORMAT}\n" \
  "    This mandatory argument specifies the format of input" \
  " frame files.  " image_file_format_arg_INFO "\n" \
  "\n" \
  "  -image {Z_k} {IMAGE_k} \n" \
  "    This argument specifies the filename of a frame from with" \
  " focus midplane at height {Z_k}.  It should be repeated for" \
  " each frame in the stack.\n" \
  "\n" \
  "  -gammaDec {GAMMADEC}\n" \
  "  -bias {BIAS}\n" \
  "    These optional arguments specify the gamma exponent and the" \
  " bias parameter to use when decoding quantized samples from the" \
  " image files into light intensity levels.  If omitted, the" \
  " values specified or implied in the input file itself are used.\n" \
  "\n" \
  "  -noise {NOISE}\n" \
  "    This argument specifies the sample noise level that one should expect in a totally out of focus image.  If not specified, it defaults to " stringify(multifok_DEFAULT_NOISE) "\n" \
  "\n" \
  "  -spread {SPREAD}\n" \
  "    This argument specifies the blur radius for spreading the pixels of the focus map.  If not specified, it defaults to " stringify(multifok_DEFAULT_SPREAD) "\n" \
  "\n" \
  "  -outPrefix {OPREF}\n" \
  "    This mandatory argument specifies the common prefix of all output files.\n" \
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
  "  Created Jun/2017 by J.Stolfi.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " multifok_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <vec.h>
#include <jsfile.h>
#include <argparser.h>
#include <wt_table.h>
#include <wt_table_gaussian.h>

#include <float_image.h>

#include <float_image_read_gen.h>
#include <float_image_write_gen.h>
#include <sample_conv.h>
#include <image_file_format.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>
#include <multifok_score.h>

#define multifok_DEFAULT_NOISE 0.001

#define multifok_DEFAULT_SPREAD 3.0

typedef struct options_t 
  { image_file_format_t imageFormat; /* Format of frame files. */
    float_vec_t zp;                  /* Z-coordinates of each image's midplane. */
    string_vec_t iname;              /* Filenames of the input images of the stack. */
    double gammaDec;                 /* Gamma exponent to assume for sample decoding, or {NAN} if omitted. */
    double bias;                     /* Bias to assume for sample decoding, or {NAN} if omitted. */
    double noise;                    /* Noise level in fully defocused image. */
    double spread;                   /* Focus spread radius. */
    char *outPrefix;                 /* Prefix for output file names. */
    bool_t verbose;                  /* TRUE to print debugging info. */
  } options_t;

options_t *multifok_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

float_image_t **multifok_read_images
  ( string_vec_t *iname, 
    image_file_format_t ffmt, 
    double gammaDec, 
    double bias, 
    bool_t verbose
  );
  /* Reads the image files whose names are {iname.e[0..n-1]} where {n =
    image.ne}, converting them from the format {ffmt} to the in-memory
    {float_image_t} format. 
    
    Applies gamma decoding with exponent {gammaDec} and specified
    {bias}. If {gammaDec} or {bias} is {NAN}, uses the exponent and bias
    specified or implied in the file. Allocates and returns a vector of
    {n} pointers to descriptors of those images. */

float_image_t *multifok_compute_focus_map
  ( float_image_t *img, 
    int32_t NW, 
    float *phi[], 
    double w[], 
    double noise,
    bool_t verbose
  );
  /* Returns a single-channel image {foc}, with same column and line
    count as {img}, such that the sample value at each pixel of {foc} is
    a measure of the sharpness of focus in the immediate neighborhood of
    the corresponding pixel of {img}.
    
    Assumes that a completely defocused image would have score {noise}. */
    
void multifok_blur_focus_map(float_image_t *foc, double spread, bool_t verbose);
  /* If {spread > 0}, the focus indicator map {foc} is blurred with an
    almost-Gaussian blur filter of deviation {spread} in each
    coordinate. */

void multifok_merge_images
  ( int32_t NI, 
    float_image_t *img[], 
    float zp[],
    float_image_t *foc[], 
    float_image_t *mrg, 
    float_image_t *zmp
  ); 
  /* Combines the images {img[0..NI-1]} into the image {mrg}, using the
    focus indicator images {foc[0..NI-1]}. Also sets channel 0 of {zmp}
    to the estimted height of each pixel, and channel 1 to the
    uncertainty in that estimate. Assumes that the focus plane of image
    {img[k]} is at height {zp[k]}.  The images {foc[k]} are 
    normalized so that their sum at every pixel is 1. */

double mutifok_focus_avg(int n, float zp[], float wt[]);
  /* Computes the average of the values {zp[0..n-1]} weighted by {wt[0.n-1]}. */
  
double mutifok_focus_var(int n, float zp[], float wt[], double avg);
  /* Computes the variance of the values {zp[0..n-1]} weighted by {wt[0.n-1]}, given their 
    nominal average. */
    
void multifok_write_merged_image
  ( float_image_t *fimg, 
    char *prefix, 
    bool_t verbose
  );
  /* Writes to "{prefix}-final.png" the image {fimg}, as a
    16-bit-per-sample PNG file with the same number of channels (1 or 3).
    Assumes that the input sample range is {[0__1]}.
    
    Applies the BT709 gamma encoding transformation to each sample. Then
    maps and quantizes the resulting samples linearly from {[0__1]} to
    {0..65535}. */

void multifok_write_all_mask_images
  ( int NI, 
    float_image_t *foc[], 
    double spr, 
    char *prefix, 
    char *tag, 
    bool_t verbose
  );
  /* Writes images {foc[0..NI-1]} to files "{prefix}-{KKK}-{tag}.png", where
    {KK} is the image index as three decimal digits, zero-padded. 
    
    Each image is written using {multifok_write_mask_image}. The
    parameter {vmax} given to {multifok_write_mask_image} will be the
    average of all samples of all images, plus {spr} times the standard
    deviation of those samples. */
    
void multifok_write_mask_image
  ( float_image_t *fimg, 
    double vmax, 
    char *prefix, 
    int32_t k, 
    char *tag, 
    bool_t verbose
  );
  /* Writes to "{prefix}-{KKK}-{tag}.png" the image {fimg} in PNG
    format; where {KKK} is the value of {k} as zero-padded 3 digits. 
    
    All samples will be scaled linearly from {[0__vmax]} to {[0__1]}. If
    {vmax} is {NAN}, assumes {vmax} equal to the max sample in {fimg}.
    Then the samples will be encoded with the BT709 gamma encoding,
    finally scaled and quantized linearly from {[0__1]} to {0..65535}. */

void multifok_write_height_map
  ( float_image_t *fimg, 
    float zmin, 
    float zmax, 
    char *prefix, 
    bool_t verbose
  );
  /* Writes to "{prefix}-height.fni" the three-channel height map {fimg}.
    
    Also writes it as a PNG file "{prefix}-height.png".
    In the latter, all samples will be scaled linearly from {[zmin__zmax]} to {[0__1]}.
    Then the samples will be encoded with the BT709 gamma encoding,
    finally scaled and quantized linearly from {[0__1]} to {0..65535}.  */

int main(int argc, char **argv)
  {
    options_t *o = multifok_parse_options(argc, argv);
    
    /* Compute the focus operator tables: */
    if (o->verbose) { fprintf(stderr, "computing focus operator tables...\n"); }
    int32_t NW = 3; /* Focus detector's window size. */
    double *w = multifok_focus_op_prod_weights(NW);
    float **phi = multifok_focus_op_basis(NW);
    /* multifok_focus_op_check(NW); */
    
    /* Get the midplane heights and number of images: */
    int NI = o->zp.ne;
    assert(NI > 0);
    float *zp = o->zp.e;
    double dzp = ((double)(zp[NI-1] - zp[0]))/NI; /* Average spacing between focus planes. */
    float zlo = (float)(zp[0] - 3.0*dzp);         /* Min height for PNG encoding. */
    float zhi = (float)(zp[NI-1] + 3.0*dzp);      /* Max height for PNG encoding. */

    /* Get the images: */
    if (o->verbose) { fprintf(stderr, "reading input images...\n"); }
    assert(o->iname.ne == NI);
    float_image_t **img = multifok_read_images
      ( &(o->iname), o->imageFormat, o->gammaDec, o->bias, o->verbose );

    int32_t NC = (int32_t)img[0]->sz[0];
    int32_t NX = (int32_t)img[0]->sz[1];
    int32_t NY = (int32_t)img[0]->sz[2];
    
    /* Compute the raw focus score images: */
    if (o->verbose) { fprintf(stderr, "computing focus score images...\n"); }
    float_image_t **foc = notnull(malloc(NI*sizeof(float_image_t*)), "no mem");
    for (int32_t k = 0; k < NI; k++)
      { foc[k] = multifok_compute_focus_map(img[k], NW, phi, w, o->noise, o->verbose); }
      
    /* Write the raw focus score images: */
    double focSpr = 2.0; /* On PNG encoding, assume max is this many standard devs above mean. */
    multifok_write_all_mask_images(NI, foc, focSpr, o->outPrefix, "foc", o->verbose);
    
    /* Spread focus information laterally: */
    if (o->verbose) { fprintf(stderr, "spreading focus scores...\n"); }
    for (int32_t k = 0; k < NI; k++)
      { multifok_blur_focus_map(foc[k], o->spread, o->verbose); }

    /* Merge images and compute the height map: */
    if (o->verbose) { fprintf(stderr, "merging images...\n"); }
    float_image_t *mrg = float_image_new(NC, NX, NY); /* Merged image. */
    float_image_t *zmp = float_image_new(2, NX, NY);  /* Height map. */
    multifok_merge_images(NI, img, zp, foc, mrg, zmp);
    
    multifok_write_merged_image(mrg, o->outPrefix, o->verbose);
    multifok_write_height_map(zmp, zlo, zhi, o->outPrefix, o->verbose);
   
    /* Write the mixing coefficient images: */
    double mixSpr = 2.0; /* On PNG encoding, assume max is this many standard devs above mean. */
    multifok_write_all_mask_images(NI, foc, mixSpr, o->outPrefix, "mix", o->verbose);
    
    if (o->verbose) { fprintf(stderr, "done.\n"); }
   
    return 0;
  }
    
void multifok_merge_images
  ( int32_t NI, 
    float_image_t *img[], 
    float zp[],
    float_image_t *foc[], 
    float_image_t *mrg, 
    float_image_t *zmp
  )
  {
    int32_t NC = (int32_t)img[0]->sz[0];
    int32_t NX = (int32_t)img[0]->sz[1];
    int32_t NY = (int32_t)img[0]->sz[2];
    
    /* Check image sizes: */
    for (int k = 0; k < NI; k++) 
      { float_image_check_size(img[k], NC, NX, NY); 
        float_image_check_size(foc[k], 1, NX, NY); 
      }
    float_image_check_size(mrg, NC, NX, NY);
    float_image_check_size(zmp, 2, NX, NY);
        
    /* Data for one pixel: */
    float *fxy = notnull(malloc(NI*sizeof(float)), "no mem"); /* Focus scores of 1 pixel. */
    float pix[NC]; /* Input image pixel. */
    float pox[NC]; /* Output image pixel. */
    for (int32_t iy = 0; iy < NY; iy++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { /* Extract focus scores for this pixel, compute total score: */
            double sum = 1.0e-200; /* Fudge factor in case all scores are zero. */
            for (int32_t k = 0; k < NI; k++)
              { fxy[k] = float_image_get_sample(foc[k], 0, ix, iy);
                assert(! isnan(fxy[k]));
                assert(fxy[k] >= 0.0);
                sum += (double)(fxy[k]);
              }

            /* Convert {fxy[0..NI-1]} to partition of unity: */
            for (int32_t k = 0; k < NI; k++) 
              { fxy[k] = (float)(((double)fxy[k])/sum); }

            /* Store back into the images: */
            for (int32_t k = 0; k < NI; k++) 
              { float_image_set_sample(foc[k], 0, ix, iy, fxy[k]); }
            
            /* Merge the images with those weights: */
            for (int32_t ic = 0; ic < NC; ic++) { pox[ic] = 0.0; }
            for (int32_t k = 0; k < NI; k++)
              { float_image_get_pixel(img[k], ix, iy, pix);
                for (int32_t ic = 0; ic < NC; ic++) 
                  { pox[ic] += (float)(fxy[k]*((double)pix[ic])); }
              }
            float_image_set_pixel(mrg, ix, iy, pox);
            
            /* Compute the mean Z and deviation of pixel: */
            double avg = mutifok_focus_avg(NI, zp, fxy);
            double var = mutifok_focus_var(NI, zp, fxy, avg);
            float_image_set_sample(zmp, 0, ix, iy, (float)avg);
            float_image_set_sample(zmp, 1, ix, iy, (float)(sqrt(var)));
          }
      }
    free(fxy);
  }

float_image_t *multifok_compute_focus_map
  ( float_image_t *img, 
    int32_t NW, 
    float *phi[], 
    double w[], 
    double noise,
    bool_t verbose
  )
  { 
    /* Get image dimensions: */
    int32_t NC = (int32_t)img->sz[0];
    int32_t NX = (int32_t)img->sz[1];
    int32_t NY = (int32_t)img->sz[2];
    
    /* Allocate output image: */
    float_image_t *foc = float_image_new(1, NX, NY);
    
    /* Compute focus scores: */
    /* Data for one sample: */
    float va[NW*NW]; /* Sample values in window for one channel. */
    for (int32_t iy = 0; iy < NY; iy++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { double sumsc = 0.0;
            for (int32_t ic = 0; ic < NC; ic++)
              { bool_t rep = TRUE; /* Set out-of-bounds samples to rearest sample. */
                float_image_get_window_samples(img, ic, ix, iy, NW, NW, rep, va);
                double fxy = multifok_focus_op_score_simple(NW, va, phi, w, noise);
                assert(! isnan(fxy));
                sumsc += fxy;
              }
            float_image_set_sample(foc, 0, ix, iy, (float)sumsc);
          }
      }

    return foc;
  }
      
void multifok_blur_focus_map(float_image_t *foc, double spread, bool_t verbose)
  { 
    if (spread > 0)
      { 
        int32_t NC = (int32_t)foc->sz[0];
        int32_t NX = (int32_t)foc->sz[1];
        int32_t NY = (int32_t)foc->sz[2];
        assert(NC == 1);
      
        double_vec_t wt = wt_table_gaussian_make(0, spread, 0.0005);
        wt_table_normalize_sum(wt.ne, wt.e);
        assert((wt.ne % 2) == 1);
        int hw = (wt.ne - 1)/2;
        if (verbose) { fprintf(stderr, "blurring with sigma = %6.4f window radius = %d\n", spread, hw); }
        float_image_t *avg = float_image_new(1, NX, NY);
        float_image_local_avg_var (foc, 0, hw, wt.e, avg, 0, NULL, 0);
        float_image_assign(foc, avg);
        float_image_free(avg);
      }
  }
  
float_image_t **multifok_read_images
  ( string_vec_t *iname, 
    image_file_format_t ffmt, 
    double gammaDec, 
    double bias,
    bool_t verbose
  )
  {
    if (verbose) { fprintf(stderr, "--- enter multifok_read_images ---\n"); }
    int32_t NI = iname->ne;
    float_image_t **img = notnull(malloc(NI*sizeof(float_image_t*)), "no mem");
    
    for (int32_t k = 0; k < NI; k++)
      { char *fname = iname->e[k];
        float v0 = 0.0;
        float vM = 1.0;
        double gammaDecFile, biasFile; /* Enconding gammaDec and bias specified or implied by file. */
        float_image_t *fimg = float_image_read_gen_named
          ( fname, ffmt, v0, vM, NULL, &gammaDecFile, &biasFile, verbose );
        int32_t NC = (int32_t)fimg->sz[0]; /* Num channels. */
        int32_t NX = (int32_t)fimg->sz[1]; /* Num columns. */
        int32_t NY = (int32_t)fimg->sz[2]; /* Num rows. */
        if (verbose) { fprintf(stderr, "read %s  NC = %d NX = %d NY = %d ...\n", fname, NC, NX, NY); }
        if (verbose) { fprintf(stderr, "gammaDec = %7.5f bias = %7.5f\n", gammaDecFile, biasFile); }
        demand((NC == 1) || (NC == 3), "images should not have alpha channels");
        img[k] = fimg;
      }
    if (verbose) { fprintf(stderr, "--- exit multifok_read_images ---\n"); }
    return img;
  }

double mutifok_focus_avg(int n, float zp[], float wt[])
  { double sum_w = 0;   /* Sum of {wt[i]}. */
    double sum_wz = 0;  /* Sum of {i*wt[i]}. */
    int i;
    for (i = 0; i < n; i++)
      { double zi = (double)(zp[i]);
        double wi = (double)(wt[i]);
        sum_w += wi;
        sum_wz += wi*zi;
      }
    double avg = sum_wz/sum_w;
    return avg;
  }
   
double mutifok_focus_var(int n, float zp[], float wt[], double avg)
  { double sum_w = 0;   /* Sum of {wt[i]}. */
    double sum_wd2 = 0;
    int i;
    for (i = 0; i < n; i++)
      { double zi = (double)(zp[i]);
        double wi = (double)(wt[i]);
        double di = zi - avg;
        sum_w += wi;
        sum_wd2 += wi*di*di; 
      }
    double var = sum_wd2/sum_w;
    return var;
  }
    
void multifok_write_merged_image
  ( float_image_t *fimg, 
    char *prefix, 
    bool_t verbose
  )
  { 
    char *fname = NULL;
    asprintf(&fname, "%s-final.png", prefix);
    double vmax = 1.0;
    if (verbose) { fprintf(stderr, "writing %s ...\n", fname); }
    float v0 = 0.0;
    float vM = 1.0;
    image_file_format_t ffmt = image_file_format_PNG;
    double gammaEnc  = sample_conv_BT709_ENC_GAMMA;
    double bias = sample_conv_BT709_BIAS;
    float_image_write_gen_named(fname, fimg, ffmt, v0, vM, gammaEnc, bias, vmax);
    free(fname);
  }

void multifok_write_all_mask_images
  ( int NI, 
    float_image_t *foc[], 
    double spr, 
    char *prefix, 
    char *tag, 
    bool_t verbose
  )
  {
    /* Choose a common scale {vmax} for the mixing coefs images: */
    if (verbose) { fprintf(stderr, "choosing scale for mixing coefs images...\n"); }
    double vmax = 1.0e-100;
    for (int32_t k = 0; k < NI; k++)
      { double avg, dev;
        float_image_compute_sample_avg_dev(foc[k], 0, &avg, &dev);
        double tmax = avg + spr*dev;
        if (tmax > vmax) { vmax = tmax; }
      }

    if (verbose) { fprintf(stderr, "writing the '%s' images (vmax = %11.8f)...\n", tag, vmax); }
    for (int32_t k = 0; k < NI; k++)
      { multifok_write_mask_image(foc[k], vmax, prefix, k, tag, verbose); }
  }

void multifok_write_mask_image
  ( float_image_t *fimg, 
    double vmax, 
    char *prefix, 
    int32_t k, 
    char *tag, 
    bool_t verbose
  )
  { 
    int32_t NC = (int32_t)fimg->sz[0];
    char *fname = NULL;

    asprintf(&fname, "%s-%03d-%s.png", prefix, k, tag);
    if (isnan(vmax))
      { /* Adjust range for this image only: */
        float sMin = 0.0f, sMax = 1.0e-35f; /* Sample range (to be computed). */
        for (int ic = 0; ic < NC; ic++) { float_image_update_sample_range(fimg, ic, &sMin, &sMax); }
        vmax = ( sMin >= 0.0 ? sMax : fmax(-(double)sMin, (double)sMax) );
      }
    if (verbose) { fprintf(stderr, "writing %s ...\n", fname); }
    float v0 = 0.0;
    float vM = (float)vmax;
    image_file_format_t ffmt = image_file_format_PNG;
    double gammaEnc = 1.0; /* Use linear encoding. */
    double bias = 0.0;
    float_image_write_gen_named(fname, fimg, ffmt, v0, vM, gammaEnc, bias, verbose);
    free(fname);
  }

void multifok_write_height_map
  ( float_image_t *fimg, 
    float zmin, 
    float zmax, 
    char *prefix, 
    bool_t verbose
  )
  { 
    /* Write the FNI version: */
    char *fname = NULL;
    asprintf(&fname, "%s-height.fni", prefix);
    if (verbose) { fprintf(stderr, "writing %s ...\n", fname); }
    FILE *wr = open_write(fname, verbose);
    float_image_write(wr, fimg);
    fclose(wr);
    free(fname);

    /* Write the PNG version: */
    asprintf(&fname, "%s-height.png", prefix);
    if (verbose) { fprintf(stderr, "writing %s ...\n", fname); }
    float v0 = zmin;
    float vM = zmax;
    image_file_format_t ffmt = image_file_format_PNG;
    double gammaEnc = 1.0; /* Use linear encoding. */
    double bias = 0.0;
    float_image_write_gen_named(fname, fimg, ffmt, v0, vM, gammaEnc, bias, verbose);
    free(fname);
 }

options_t *multifok_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    /* Parse keyword parameters: */

    o->imageFormat = image_file_format_arg_parse(pp, "-imageFormat");

    if (argparser_keyword_present(pp, "-gammaDec"))
      { o->gammaDec = argparser_get_next_double(pp, 0.2, 5.0); }
    else
      { o->gammaDec = NAN; }

    if (argparser_keyword_present(pp, "-bias"))
      { o->bias = argparser_get_next_double(pp, 0.0, 0.5); }
    else
      { o->bias = NAN; }

    if (argparser_keyword_present(pp, "-noise"))
      { o->noise = argparser_get_next_double(pp, 0.000, 1.000); }
    else
      { o->noise = multifok_DEFAULT_NOISE; }

    if (argparser_keyword_present(pp, "-spread"))
      { o->spread = argparser_get_next_double(pp, 0.000, 100.000); }
    else
      { o->spread = multifok_DEFAULT_SPREAD; }

    o->zp = float_vec_new(100);
    o->iname = string_vec_new(100);
    int32_t NI = 0;
    while (argparser_keyword_present(pp, "-image"))
      { float_vec_expand(&(o->zp), NI);
        o->zp.e[NI] = (float)argparser_get_next_double(pp, -1000.0, +1000.0);
        string_vec_expand(&(o->iname), NI);
        o->iname.e[NI] = argparser_get_next_non_keyword(pp); 
        NI++;
      }
    if (NI == 0) { argparser_error(pp, "must specifiy at least one \"-image\""); }
    float_vec_trim(&(o->zp), NI);
    string_vec_trim(&(o->iname), NI);
        
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next(pp);

    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }


