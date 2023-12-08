#define PROG_NAME "multifok_make_stack"
#define PROG_DESC "Creates a synthetic multi-focus image stack"
#define PROG_VERS "1.0"

// Last edited on 2023-11-26 06:38:29 by stolfi

#define multifok_make_stack_C_COPYRIGHT \
    "© 2018 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -imageFormat {IMAGE_FORMAT} \\\n" \
  "    -image {TIMG} \\\n" \
  "    -heightMap {ZMAP} \\\n" \
  "    -numFrames {NZ} \\\n" \
  "    -minHeight {ZMIN} \\\n" \
  "    -heightStep {ZSTEP} \\\n" \
  "    -blurFactor {BLF} \\\n" \
  "    [ -minBlur {BLMIN} ] \\\n" \
  "    -outDir {OUT_DIR} \\\n" \
  "    -framePattern {FRAME_PATTERN} \\\n" \
  "    [ -verbose ] \\\n"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "   Reads an image {TIMG} that is sharp everywhere, and a height map {ZMAP} that specifies the {Z}-coordinate of the object's visible surface in each pixel.  Then computes and writes out a stack of images (/frames/) {FIMG[0..NZ-1]} that simulate the result of photographing the image with a camera that has a" \
  " limited depth of focus, different object-camera" \
  " distances.\n" \
  "\n" \
  "   The output frames will be numbered in {0..NZ-1}.  The {Z}-coordinate {ZF[k]} of the focus plane of frame {FIMG[k]} is {ZMIN + k*ZSTEP}.\n" \
  "\n" \
  "   The input files must be in the specified format, both with" \
  " the same size, and must be aligned.  On input, each pixel of the height map {ZMAP} is linearly mapped from its range {0..MAXVL} in the file to a real number in {[0 _ 1]}." \
  "\n" \
  "   Each pixel {TIMG[p]} of {TIMG} is turned into a Gaussian blurred spot, centered on puxel {p}, and added to each frame {FIMG[k]}.  The mean radius of the blurred spot is approximately {{BLF]*|ZMAP[p] - ZF[k]|}.  However, this formula is modified so that the blur radius is never less than {BLMIN}, even when {ZMAP[p]} is equal to {ZF[k]}.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  For each integer {k} in {0..NZ-1}, the program writes an image" \
  " file \"{OUT_DIR}/{KKKKK}.png\", where {KKKKK} is the value of {k} formatted with {FRAME_PATTERN}.\n" \
  "\n" \
  "  For each output image above, the program also writes a mask image" \
  " file \"{OUT_DIR}/{KKKKK}_mask.png\", which (if converted to float) is 1.0 where the frame is perfectly in focus, and 0.0 where it is aximally unfocused.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -imageFormat {IMAGE_FORMAT}\n" \
  "    This mandatory argument specifies the format of input" \
  " image files.  " image_file_format_arg_INFO "\n" \
  "\n" \
  "  -image {TIMG} \n" \
  "    This mandatory argument specifies the file name of the input image, with" \
  " the object entirely in focus.\n" \
  "\n" \
  "  -heightMap {ZMAP} \n" \
  "    This mandatory argument specifies the file name of the input height map.\n" \
  "\n" \
  "  -numFrames {NZ} \n" \
  "  -minHeight {ZMIN} \n" \
  "  -heightStep {ZSTEP} \n" \
  "    These mandatory arguments specify the number {NZ} of frames to generate, the {Z}-coordinate {ZMIN} of the focus midplane of frame 0, and the {Z} displacement {ZSTEP} between successive frames.\n" \
  "\n" \
  "  -blurFactor {BLF} \n" \
  "    This mandatory argument specifies the proportionality factor between the {Z}-distance between the surface and the frame's focus midplaneand blurring radius.\n" \
  "\n" \
  "  -minBlur {BLMIN} \n" \
  "    This optional argument specify the minimum blurring radius, when the surface is at the focus midplane.  If omitted, the program assumes {BLMIN=0}.\n" \
  "\n" \
  "  -outDir {OUT_DIR}\n" \
  "    This mandatory argument specifies the directory were all output files will be written to.\n" \
  "\n" \
  "  -framePattern {FRAME_PATTERN}\n" \
  "    This mandatory argument specifies the pattern for frame file names.  It must include exacly one unquoted {printf} formatting code for a 32-bit signed integer, like '%06d'.  That code will be replaced by the frame number {k}.  The {FRAME_PATTERN} must NOT include the file extension.\n" \
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
  "  Created Mar/2018 by J.Stolfi.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2018-03-02 Created. [J.Stolfi]\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " multifok_make_stack_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#define _GNU_SOURCE
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
#include <wt_table.h>

#include <float_image.h>
#include <float_image_paint.h>
#include <float_image_read_gen.h>
#include <float_image_write_gen.h>
#include <sample_conv.h>
#include <image_file_format.h>

#include <multifok_make_stack_read_image.h>

typedef struct options_t 
  { image_file_format_t imageFormat; /* Format of frame files. */
    char *image;                     /* Filename of the perfectly focused image. */
    char *heightMap;                 /* Filename of the height map. */
    int32_t numFrames;               /* Number of frames to generate. */
    double minHeight;                /* {Z}-coordinate of frame 0. */
    double heightStep;               /* {Z} increment between frames. */
    double blurFactor;               /* Factor in blur radius formula. */
    double minBlur;                  /* Min blur radius. */
    char *outDir;                    /* Directory for output files. */
    char *framePattern;              /* Filename pattern of the output frames.. */
    bool_t verbose;                  /* TRUE to print debugging info. */
  } options_t;

options_t *multifok_make_stack_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */

float_image_t **multifok_make_stack_alloc_images(int32_t NF, int32_t NX, int32_t NY);
  /* ALlocates {NF} monchrome images with {NX} columns and {NY} rows of pixels,
    and fills them with zeros. */

void multifok_make_stack_compute_blurred_images
  ( float_image_t *timg, 
    float_image_t *zimg, 
    int32_t NF, 
    double minHeight,                /* {Z}-coordinate of frame 0. */
    double heightStep,               /* {Z} increment between frames. */
    double blurFactor,               /* Factor in blur radius formula. */
    double minBlur,                  /* Min blur radius. */
    bool_t verbose, 
    float_image_t **frame,           /* (OUT) frames. */
    float_image_t **fMask            /* (OUT) focus masks. */
  );
  /* Fills {frame[0..NF-1]} with {NF} blurred versions of {timg}, corresponding to camera
    positions that place the focus plane at {Z}-coordinates described by
    {minHeight} and {heightStep}. The blurring of each pixel of {timg}
    is defined by its {Z}-coordinate in {zimg}, and parameters
    {blurFactor} and {minBlur}.  Also stores in {fMask[0..NF-1]} 
    the focus mask for each image. 
    
    The images {frame[0..NF-1]} and *{fMask[0..NF-1]} must be allocated 
    by the caller. */

void multifok_make_stack_compute_blurred_image
  ( float_image_t *timg, 
    float_image_t *zimg, 
    double frameZ,                /* {Z}-coordinate of frame. */
    double blurFactor,            /* Factor in blur radius formula. */
    double minBlur,               /* Min blur radius. */
    float_image_t *frame,         /* (OUT) frame. */
    float_image_t *fMask          /* (OUT) focus mask. */
  );
  /* Stores in {frame} a blurred version of {timg}, corresponding to camera
    placement that puts the focus plane at {Z}-coordinate {frameZ}. The
    blurring of each pixel of {timg} is defined by its {Z}-coordinate in
    {zimg}, and parameters {blurFactor} and {minBlur}. Also stores in {fMask} 
    the focus mask. */

void multifok_make_stack_write_frame_images
  ( char *outDir, 
    int32_t NF, 
    char *framePattern,              /* Filename pattern of the output frames. */
    bool_t verbose, 
    float_image_t **frame, 
    float_image_t **fMask
  );
  /* Assumes that {NF} is the number of input frames, and that {frame}
    is a vector of {NF} pointers to images, all monochromatic with the
    same size as the input frames, with samples in {[0_1]}.
    
    Writes each image {frame[f]}  to files "{outDir}/frame_{NNNNN}.png" where {NNNNN}
    is the numeric frame ID {frameID[f]} formatted as 5 digits, zero-padded. */

void multifok_make_stack_analyze_image(float_image_t *img);
  /* Report the statistics of image {img} to {stderr}. */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = multifok_make_stack_parse_options(argc, argv);
    
    /* Get the input image: */
    if (o->verbose) { fprintf(stderr, "reading input image...\n"); }
    float_image_t *timg = multifok_make_stack_read_image(o->image, o->imageFormat, o->verbose);
    int32_t NC, NX, NY;
    float_image_get_size(timg, &NC, &NX, &NY);
    demand(NC == 1, "input image must be monochromatic");
    
    /* Get the height map: */
    if (o->verbose) { fprintf(stderr, "reading height map...\n"); }
    float_image_t *zimg = multifok_make_stack_read_image(o->heightMap, o->imageFormat, o->verbose);
    float_image_check_size(zimg, 1, NX, NY);

    /* Allocate the output frames and masks, fill with zeros: */
    if (o->verbose) { fprintf(stderr, "allocating the output images...\n"); }
    int32_t NF = o->numFrames;
    float_image_t **frame = multifok_make_stack_alloc_images(NF, NX, NY);
    float_image_t **fMask = multifok_make_stack_alloc_images(NF, NX, NY);
    
    /* Compute the blurred images: */
    if (o->verbose) { fprintf(stderr, "computing the blurred images and masks...\n"); }
    multifok_make_stack_compute_blurred_images
      ( timg, zimg, 
        NF, o->minHeight, o->heightStep, o->blurFactor, o->minBlur, 
        o->verbose, 
        frame, fMask
      );
    
    /* Write the similarity images: */
    multifok_make_stack_write_frame_images(o->outDir, NF, o->framePattern, o->verbose, frame, fMask);

    if (o->verbose) { fprintf(stderr, "done.\n"); }
    return 0;
  }

float_image_t **multifok_make_stack_alloc_images(int32_t NF, int32_t NX, int32_t NY)
  { float_image_t **img = notnull(malloc(NF*sizeof(float_image_t *)), "no mem");
    for (int32_t i = 0; i < NF; i++)
      { img[i] = float_image_new(1, NX, NY);
        float_image_fill_channel(img[i], 0, 0.0f);
      }
    return img;
  }

void multifok_make_stack_compute_blurred_images
  ( float_image_t *timg, 
    float_image_t *zimg, 
    int32_t NF, 
    double minHeight,                /* {Z}-coordinate of frame 0. */
    double heightStep,               /* {Z} increment between frames. */
    double blurFactor,               /* Factor in blur radius formula. */
    double minBlur,                  /* Min blur radius. */
    bool_t verbose, 
    float_image_t **frame,           /* (OUT) frames. */
    float_image_t **fMask            /* (OUT) focus masks. */
  )
  {
    for(int32_t f = 0; f < NF; f++)
      { double frameZ = minHeight + f*heightStep;
        if (verbose) { fprintf(stderr, "computing frame %d (Z = %.5f)...\n", f, frameZ); }
        multifok_make_stack_compute_blurred_image(timg, zimg, frameZ, blurFactor, minBlur, frame[f], fMask[f]);
      }
   }

void multifok_make_stack_compute_blurred_image
  ( float_image_t *timg, 
    float_image_t *zimg, 
    double frameZ,                /* {Z}-coordinate of frame. */
    double blurFactor,            /* Factor in blur radius formula. */
    double minBlur,               /* Min blur radius. */
    float_image_t *frame,         /* (OUT) frame. */
    float_image_t *fMask          /* (OUT) focus mask. */
  )
  { 
    int32_t NC, NX, NY;
    float_image_get_size(timg, &NC, &NX, &NY);
    /* Paranoia: */
    demand(NC == 1, "input image must be monochromatic");
    float_image_check_size(zimg, 1, NX, NY);
    
    /* Clear the images: */
    float_image_fill(frame, 0.0f);
    float_image_fill(fMask, 0.0f);
    
    int32_t subm = 1; /* Subsampling order. */
    bool_t debugged = FALSE;
    for (int32_t iy = 0; iy < NY; iy++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { /* Pixel center coordinates: */
            double xctr = (double)ix + 0.5;
            double yctr = (double)iy + 0.5;
            /* Pixel value in the original image: */
            double tv = (double)float_image_get_sample(timg, 0, ix, iy);
            /* Compute the blur radius {rad}: */
            float fz = float_image_get_sample(zimg, 0, ix, iy);
            double dz = ((double)fz) - frameZ;
            double rad = hypot(blurFactor*dz, minBlur);
            /* Compute the integral of the smudge for unit intensity: */
            double vTot = float_image_paint_smudge(NULL, 0, xctr, yctr, rad, rad, 1.0f, subm);
            /* Splat a smudge with total mass {tv}: */
            float fv = (float)(tv / vTot);
            if ((! debugged) && (fv > 0.1))
              { fprintf(stderr, "fz = %.5f dz = %.5f rad = %4f", fz, dz, rad);
                fprintf(stderr, " vTot = %.4f tv = %.4f fv = %.4f\n", vTot, tv, fv);
                debugged = TRUE;
              }
            (void)float_image_paint_smudge(frame, 0, xctr, yctr, rad, rad, fv, subm);
            /* Set the mask image: */
            float fm = (float)(1.0/hypot(1.0, rad));
            float_image_set_sample(fMask, 0, ix, iy, fm);
          }
      }
  }

void multifok_make_stack_write_frame_images
  ( char *outDir, 
    int32_t NF, 
    char *framePattern,              /* Filename pattern of the output frames. */
    bool_t verbose, 
    float_image_t **frame, 
    float_image_t **fMask
  )
  {
    for(int32_t f = 0; f < NF; f++)
      { char *name = NULL;
        asprintf(&name, "%s/frame_%05d", outDir, f);
        image_file_format_t ffmt = image_file_format_PNG;
        double v0 = 0.0;
        double vM = 1.0;
        double gammaEnc = 1.0;
        double bias = 0.0;
        bool_t verbose = FALSE;
        char *fname = NULL;
        
        /* Write the frame:  */
        asprintf(&fname, "%s.png", name);
        fprintf(stderr, "writing %s\n", fname);
        fprintf(stderr, "channel ranges before gamma encoding:\n");
        multifok_make_stack_analyze_image(frame[f]);
        float_image_write_gen_named(fname, frame[f], ffmt, v0, vM, gammaEnc, bias, verbose);
        free(fname);
        
        /* Write the mask: */
        asprintf(&fname, "%s_mask.png", name);
        fprintf(stderr, "writing %s\n", fname);
        float_image_write_gen_named(fname, fMask[f], ffmt, v0, vM, gammaEnc, bias, verbose);
        free(fname);
        free(name);
      }
  }
    
void multifok_make_stack_analyze_image(float_image_t *img)
  { int32_t NC = (int32_t)img->sz[0];;
    for (int32_t c = 0; c < NC; c++)
      { float vMin = +INF;
        float vMax = -INF;
        float_image_update_sample_range(img, c, &vMin, &vMax);
        fprintf(stderr, "  channel %d pixel range = [ %.4f __ %.4f ]\n", c, vMin, vMax);
      }
  }

options_t *multifok_make_stack_parse_options(int32_t argc, char **argv)
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

    argparser_get_keyword(pp, "-image");
    o->image = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-heightMap");
    o->heightMap = argparser_get_next_non_keyword(pp);

    argparser_get_keyword(pp, "-numFrames");
    o->numFrames = (int32_t)argparser_get_next_int(pp, 1, 1000);
        
    argparser_get_keyword(pp, "-minHeight");
    o->minHeight = argparser_get_next_double(pp, -2.00, +2.00);
        
    argparser_get_keyword(pp, "-heightStep");
    o->heightStep = argparser_get_next_double(pp, 0.0001, 4.000);
        
    argparser_get_keyword(pp, "-blurFactor");
    o->blurFactor = argparser_get_next_double(pp, 0.01, 10.000);

    if (argparser_keyword_present(pp, "-minBlur"))
      { o->minBlur = argparser_get_next_double(pp, 0.000, 10.000); }
    else
      { o->minBlur = 0.0; }

    argparser_get_keyword(pp, "-outDir");
    o->outDir = argparser_get_next(pp);

    argparser_get_keyword(pp, "-framePattern");
    o->framePattern = argparser_get_next_non_keyword(pp);

    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }


