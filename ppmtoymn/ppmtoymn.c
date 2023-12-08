#define PROG_NAME "ppmtoymn"
#define PROG_DESC "convert PPM image from RGB space to an luv-type (ymn) space"
#define PROG_VERS "1.0"

/* Last edited on 2021-12-31 23:43:52 by stolfi */

#define ppmtoymn_C_COPYRIGHT "Copyright © 2003 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    [ -rgbscale {RWT} {GWT} {BWT} ] \\\n" \
  "    [ -uvrelative {UVREL} ] \\\n" \
  "    [ -yuvscale {YWT} {UWT} {VWT} ] \\\n" \
  "    [ -uvcompand {UVCOMP} ] [ -ycompand {YCOMP} ] \\\n" \
  "    [ -inverse ] \\\n" \
  "    [ PPMFILE ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads the PPMFILE, and converts the pixel values from" \
  " RGB to a different color space, called here ymn, derived" \
  " from YUV (roughly, brightness/redness/blueness), or" \
  " vice-versa.  The new image is written to the standard output.\n" \
  "  Optionally, the {U} and {V} coordinates may be divided" \
  " by {Y}, so as to express hue and saturation independently" \
  " of the brightness. Other options include arbitrary" \
  " scaling of the RGB and the YMN coordinates, and nonlinear" \
  " mapping of Y and/or U,V  to compensate for quantization error.\n" \
  "  To be precise, the transformation consists of the following steps:\n" \
  "\n" \
  "    0. Map the input RGB pixels from integers in {0..maxval} to reals in [0,1].\n" \
  "    1. Unscale (divide) the {R,G,B} coords by {RWT,GWT,BWT}.\n" \
  "    2. Linearly map the {R,G,B} coordinates to {Y,U,V}.\n" \
  "    3. Make {U,V} relative to {Y} with {UVREL} bias, to get {u,v}.\n" \
  "    4. Divide the {u,v} coordinates by {uvnorm}.\n" \
  "    5. Scale (multiply) the {Y,u,v} coordinates by {YWT,UWT,VWT}.\n" \
  "    6. Clip {u,v} to the unit disk, preserving hue; and {Y} to [0_1].\n" \
  "    7. Compand {Y} with the {YCOMP} parameter, obtaining {y}.\n" \
  "    8. Compand radially {u,v} with {UVCOMP} parameter, obtaining {m,n}.\n" \
  "    9. Linearly quantize the {y} coordinate, from the" \
  " range [0_1] to {0..maxval}; and the the {m} and {n}" \
  " coordinates, from the   range [-1_+1] to {0.. maxval}.\n" \
  "\n" \
  "  Step 1 is a linear transformation from {R,G,B} to {Y,U,V} coordinates" \
  " where {Y} is the luminance" \
  " (brightness)  in the range [0_1], and {U}, {V} are" \
  " luminance-free chroma coordinates, roughly" \
  " \"redness\" and \"blueness\", each in the range [-1_+1]." \
  "  This step uses {frgb_to_YUV_b} in {frgb_ops.h}.\n" \
  "\n" \
  "  Step 2 leaves {Y} in [0_1].\n" \
  "\n" \
  "  Step 3 consists of scaling {U,V} by {s = (1 + a)/(Y + a)}, where\n" \
  " {a = UVREL}. If {UVREL} is {+INFINITY}, {s} is 1 and step 3 is a no-op.\n" \
  "\n" \
  "  In step 5, the {Y} coordinates is scaled relative to 0.5 rather\n" \
  " than 0; that is, by doing {Y = YWT*(Y - 0.5) + 0.5}.  Thus,\n" \
  " setting {YWT=0} leaves {Y=0.5} for any color.\n" \
  "\n" \
  "  Step 7 consists of {Y = compand(Y*ycompA, YCOMP, FALSE)/ycompB}; where\n" \
  " {ycompB = compand(1.0/ycompA, YCOMP)}, so that [0_1] gets mapped to [0_1].\n" \
  " If {YCOMP} is {1.0}, step 7 is a no-op.\n" \
  "\n" \
  "  Step 8 consists of multiplying {u/r,v/r} by\n" \
  " {compand(r*uvcompA,UVCOMP,FALSE)/uvcompB} where {r = sqrt(u^2+v^2)} and\n" \
  " {uvcompB = compand(1.0/uvcompA, UVCOMP)}; so that the unit {u,v} disk is\n" \
  " mapped to itself. If {UVCOMP} is {1.0}, step 8 is a no-op.\n" \
  "\n" \
  "  The {compand} function above is a cube-root-like transformation.\n" \
  "\n" \
  "  Note that the weights {RWT,GWT,YWT} are divided into" \
  " the RGB coordinates (and hence must be nonzero);" \
  " while {YWT,UWT,VWT} are multiplied into {y,u,v}" \
  " (so they may be zero).\n" \
  "\n" \
  "  The \"-inverse\" option performs the inverse" \
  " mapping, from ymn to RGB. Steps 1--8 are  inverted" \
  " and executed in the opposite order.  Thus the {YWT,UWT,VWT}" \
  " weights must be nonzero when computing the inverse.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -blabber {AMOUNT}\n" \
  "    Blabbers for that {AMOUNT}. May also bla" \
  " bla bla bla bla bla bla bla bla bla bla bla bla" \
  " bla bla bla bla bla bla bla.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ppm(5)\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ppmtoyuv(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in feb/2003 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2006-11-22 General revamp by J. Stolfi, IC-UNICAMP.\n" \
  "  2007-10-26 Updated arg parsing, logic by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " ppmtoymn_C_COPYRIGHT "\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

/* We need to set these in order to get {asinh}. What a crock... */
#undef __STRICT_ANSI__
#define _ISOC99_SOURCE 1
#include <math.h>

#include <jspnm.h>
#include <uint16_image.h>
#include <argparser.h>
#include <frgb.h>
#include <frgb_ops.h>

#include <assert.h>
#include <values.h>

typedef struct rgb_pixel_t { uint16_t c[3]; } rgb_pixel_t; 

typedef struct transform_t
  { 
    double RWT, GWT, BWT;  /* RGB input scaling weights. */
    double UVREL;       /* Bias for converting linear UV to relative uv. */
    double uvnorm;      /* Normalization factor for {u,v} */
    double YWT, UWT, VWT;  /* Ymn Output scaling weights. */
    double YCOMP;       /* Bias for Y companding. */
    double ycompA;      /* Pre-scaling factor for {Y} companding. */
    double ycompB;      /* Post-scaling factor for {Y} companding. */
    double UVCOMP;      /* Bias for {u,v} companding. */
    double uvcompA;     /* Pre-scaling factor for {u,v} companding. */
    double uvcompB;     /* Post-scaling factor for {u,v} companding. */
  } transform_t;
  /*
    Describes a colorspace transformation to be applied to 
    an {frg_t} value. */

typedef struct options_t
  { transform_t *tr;
    bool_t inverse;
    char *imgname;
  } options_t;

/* PROTOTYPES */

double compand(double v, double vcomp);
  /* Computes a nonlinear function of {v} that is monotonic, odd, and
    smooth. The parameter {vcomp} controls the derivative near the
    origin. {compand(v,1/vcomp)} gives the inverse function. */
    
void rgb_to_ymn(rgb_pixel_t *pp, uint16_t maxval, transform_t *tr);
  /* Converts {*pp} from RGB to ymn space (yuv + chroma compression)
    by the transformation parameters {tr}.
    
    On entry, the components of rgb_pixel_t {pp} are lineraly mapped
    from integers in {0..maxval} to floats in the range [0_1].
    
    On output, the {y,m,n} components are mapped linearly to the
    interval {[0 _ maxval]} and rounded to integers. Assumes that
    the range of {y} is [0_1], and that of {m,n} is [-1_+1].
  */

void ymn_to_rgb(rgb_pixel_t *pp, uint16_t maxval, transform_t *tr);
  /* Converts {*pp} from ymn space to RGB, by reversing the 
    transformation of {rgb_to_ymn}. */

double compute_uvnorm(double UVREL);
  /* Computes the maximum chroma {sqrt(u*u+v*v)} of any color
    in the RGB unit cube after transformation to {YUV}
    space and then to {Yuv} space with the given UVREL. */

void compute_compand_factors(double comp, double *compA, double *compB);
  /* Computes the (multiplicative) pre-scaling factor and the 
    (dividing) post-scale factor for the companding formula. */

void transform_ppm_image (uint16_image_t *img, transform_t *tr, bool_t inverse);
  /* Applies {rgb_to_ymn} to all pixels of {img}. */

options_t *parse_options(int argc, char **argv);
  /* Obtains arguments from the command line. */

int main(int argc, char* argv[]);

/* IMPLEMENTATIONS */

int main(int argc, char* argv[])
  { 
    /* Command-line parameters: */
    options_t *o = parse_options(argc, argv);
    
    uint16_image_t *img;

    img = uint16_image_read_pnm_named((o->imgname == NULL ? "-" : o->imgname), FALSE);
    
    transform_ppm_image(img, o->tr, o->inverse);

    uint16_image_write_pnm_named("-", img, FALSE, FALSE);
    return 0;
  }

double compute_uvnorm(double UVREL)
  {
    double rmax = 0.0;
    frgb_t p = (frgb_t){{ 0.0, 0.0, 0.0 }};
    /* The max chroma in the RGB unit cube should be at a corner: */
    int i, j;
    for (i = 0; i < 3; i++)
      { for (j = 0; j < 2; j++)
          { /* Generate a primary or secondary color: */
            p.c[0] = p.c[1] = p.c[2] = (float)j;
            p.c[i] = 1.0 - p.c[i];
            /* Conver to YUV: */
            frgb_to_YUV_b(&p);
            /* Apply non-linear step: */
            if (UVREL != +INFINITY) { frgb_YUV_to_Yuv(&p, UVREL); }
            /* Get chroma: */
            double u = p.c[1];
            double v = p.c[2];
            double rad = sqrt(u*u + v*v);
            if (rad > rmax) { rmax = rad; }
          }
      }
    return rmax;
  }

void compute_compand_factors(double comp, double *compA, double *compB)
  {
    double s = sqrt(comp);
    (*compA) = sinh(1/s);
    (*compB) = sinh(s);
  }

options_t *parse_options(int argc, char **argv)
  {
    /* Initialize the argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    /* Create a transform record, fill it with default parameters: */
    transform_t *tr = (transform_t *)notnull(malloc(sizeof(transform_t)), "no mem");

    o->imgname = NULL;
    
    if (argparser_keyword_present(pp, "-rgbscale"))
      { tr->RWT = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        tr->GWT = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        tr->BWT = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
      }
    else
      { tr->RWT = 1.0;
        tr->GWT = 1.0;
        tr->BWT = 1.0;
      }
    
    if (argparser_keyword_present(pp, "-uvrelative"))
      { tr->UVREL = argparser_get_next_double(pp, 0.0, +DBL_MAX); }
    else
      { tr->UVREL = INFINITY; }
    
    if (argparser_keyword_present(pp, "-yuvscale"))
      { tr->YWT = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
        tr->UWT = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
        tr->VWT = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); 
      }
    else
      { tr->YWT = 1.0;
        tr->UWT = 1.0;
        tr->VWT = 1.0;
      }
    
    if (argparser_keyword_present(pp, "-ycompand"))
      { tr->YCOMP = argparser_get_next_double(pp, 0.0, +DBL_MAX); }
    else
      { tr->YCOMP = 1.0; }
      
    if (argparser_keyword_present(pp, "-uvcompand"))
      { tr->UVCOMP = argparser_get_next_double(pp, 0.0, +DBL_MAX); }
    else
      { tr->UVCOMP = 1.0; }
      
    o->inverse = argparser_keyword_present(pp, "-inverse");

    /* Get positional arguments: */
    argparser_skip_parsed(pp);
    if (argparser_next(pp) != NULL)
      { o->imgname = argparser_get_next(pp); }
    else
      { o->imgname = "-"; }

    /* Check for leftover args: */
    argparser_finish(pp);

    /* Compute normalization factors: */
    tr->uvnorm = compute_uvnorm(tr->UVREL);
    compute_compand_factors(tr->UVCOMP, &(tr->uvcompA), &(tr->uvcompB));
    compute_compand_factors(tr->YCOMP, &(tr->ycompA), &(tr->ycompB));
    
    o->tr = tr;

    return o;
  }

void transform_ppm_image(uint16_image_t *img, transform_t *tr, bool_t inverse)
{ int row, col, c;
    assert(img->chns == 3);
    for (row = 0; row < img->rows; row++)
      { uint16_t *sP = img->smp[row];
        for (col = 0; col < img->cols; col++)
          { rgb_pixel_t pp;
            for (c = 0; c < img->chns; c++) { pp.c[c] = sP[c]; }
            if (inverse)
              { ymn_to_rgb(&pp, img->maxval, tr); }
            else
              { rgb_to_ymn(&pp, img->maxval, tr); }
            for (c = 0; c < img->chns; c++) { sP[c] = pp.c[c]; }
            sP += img->chns;
          }
      }
  }

void rgb_to_ymn(rgb_pixel_t *pp, uint16_t maxval, transform_t *tr)
  { /* Convert RGB to [0_1], un-applying RGB scales: */
    double m = (double)maxval;
    frgb_t p;
    p.c[0] = ((double)pp->c[0])/m/tr->RWT;
    p.c[1] = ((double)pp->c[1])/m/tr->GWT;
    p.c[2] = ((double)pp->c[2])/m/tr->BWT;
    /* Convert to yuv system: */
    frgb_to_YUV_b(&p);
    if (tr->UVREL != +INFINITY) { frgb_YUV_to_Yuv(&p, tr->UVREL); }
    /* Apply {uvnorm} and the yuv scales: */
    double wy = tr->YWT*(p.c[0] - 0.5) + 0.5;
    double wu = tr->UWT*p.c[1]/tr->uvnorm;
    double wv = tr->VWT*p.c[2]/tr->uvnorm;
    /* Clip {wy} to range [0_1]: */
    if (wy < 0.0) { wy = 0.0; } else if (wy > 1.0) { wy = 1.0; }
    /* Nonlinear companding for {wy}: */
    if ((tr->YCOMP != 1.0) && (wy != 0.0) && (fabs(wy) != 1.0))
      { wy *= tr->ycompA;
        wy = sinh(asinh(wy)*tr->YCOMP);
        wy /= tr->ycompB;
      }
    /* Clip {wu,wv} to max unit norm, preserving hue: */
    double r2 = wu*wu + wv*wv;
    if (r2 > 1.0) { double s = 1.0/sqrt(r2); wu *= s; wv *= s; r2 = 1.0; }
    /* Nonlinear companding for {wu,wv}: */
    if ((tr->UVCOMP != 1.0) && (r2 != 0.0) && (r2 != 1.0))
      { double s = r2;
        s *= tr->uvcompA;
        s = sinh(asinh(s)*tr->UVCOMP);
        s /= tr->uvcompB;
        s /= r2;
        wu *= s; wv *= s;
      }
    /* Quantization (unsigned for {wy}, signed for {wu,wv}): */
    { int iy = (int)(m*wy + 0.5);
      int iu = (int)(m*(wu + 1.0)/2.0 + 0.5);
      int iv = (int)(m*(wv + 1.0)/2.0 + 0.5);
      pp->c[0] = iy;
      pp->c[1] = iu;
      pp->c[2] = iv;
    }
  }

void ymn_to_rgb(rgb_pixel_t *pp, uint16_t maxval, transform_t *tr)
  { /* Undo quantization (unsigned for {wy}, signed for {wu,wv}): */
    double m = (double)maxval;
    double wy =       ((double)pp->c[0])/m;
    double wu = 2.0 * ((double)pp->c[1])/m - 1.0;
    double wv = 2.0 * ((double)pp->c[2])/m - 1.0;
    /* Reverse nonlinear companding of {wu,wv}: */
    double r2 = wu*wu + wv*wv;
    if ((tr->UVCOMP != 1.0) && (r2 != 0.0) && (r2 != 1.0))
      { double s = r2;
        s *= tr->uvcompB;
        s = sinh(asinh(s)/tr->UVCOMP);
        s /= tr->uvcompA;
        s /= r2;
        wu *= s; wv *= s;
      }
    /* Reverse nonlinear companding of {wy}: */
    if ((tr->YCOMP != 1.0) && (wy != 0.0) && (fabs(wy) != 1.0))
      { wy *= tr->ycompB;
        wy = sinh(asinh(wy)/tr->YCOMP);
        wy /= tr->ycompA;
      }
    /* Un-apply yuv scales: */
    frgb_t p;
    p.c[0] = (float)0.5 + (wy - 0.5)/tr->YWT;
    p.c[1] = tr->uvnorm*(float)wu/tr->UWT;
    p.c[2] = tr->uvnorm*(float)wv/tr->VWT;
    /* Convert to RGB, apply RGB scales: */
    if (tr->UVREL != +INFINITY) { frgb_YUV_from_Yuv(&p, tr->UVREL); }
    frgb_from_YUV_b(&p);
    double R = (double)p.c[0] * tr->RWT;
    double G = (double)p.c[1] * tr->GWT;
    double B = (double)p.c[2] * tr->BWT;
    /* Clip to unit cube preserving hue: */
    { double dR = R - 0.5;
      double dG = G - 0.5;
      double dB = B - 0.5;
      double dM = fabs(dR);
      if (fabs(dG) > dM) { dM = fabs(dG); }
      if (fabs(dB) > dM) { dM = fabs(dB); }
      if (dM > 0.5)
        { double s = 0.5/dM; 
          R = 0.5 + s*dR; G = 0.5 + s*dG; B = 0.5 + s*dB;
        }
    }
    if ((fabs(R) > 1.00001) || (fabs(G) > 1.00001) || (fabs(B) > 1.00001)) 
      { fprintf(stderr, "×"); }
    /* Quantization: */
    { int ir = (int)(m*R + 0.5);
      int ig = (int)(m*G + 0.5);
      int ib = (int)(m*B + 0.5);
      pp->c[0] = ir;
      pp->c[1] = ig;
      pp->c[2] = ib;
    }
  }
    
