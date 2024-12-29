#define PROG_NAME "ppmmarpap"
#define PROG_DESC "Simulate marbleized paper."
#define PROG_VERS "1.0"
/* Last edited on 2024-12-21 11:58:51 by stolfi */

#define ppmmarpap_C_COPY \
  "Copyright © 2012 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -size {NX} {NY} \\\n" \
  "  [ -source {PNM_FILE} ] \\\n" \
  "  -whorls {NW} \\\n" \
  "  [ -detailSize {NUM} ] \\\n" \
  "  [ -maxval {MAXVAL} ] \\\n" \
  "  [ -verbose ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program simulates the creation of marbleized paper, in which" \
  " a multicolor layer of oil-based paint is floated on a tray filled with water, roughly stired," \
  " then transferred to a sheet of paper.  The resulting image is writen" \
  " to the standard output\n" \
  "OPTIONS\n" \
  "  -size {NX} {NY} \n" \
  "    These mandatory parameters are the dimensions of the outout image in pixels.\n" \
  "\n" \
  "  -source {PNM_FILE} \n" \
  "    This optional parameter specifies the PGM or PPM image to be used as" \
  " the initial state of the paint layer, before the rough" \
  " stirring.  If {PNM_FILE} is \"-\", the initial image is read from the standard input. If not" \
  " specified, the program will use a series of fuzzy horizontal stripes with dark" \
  " bubble spots.\n" \
  "\n" \
  "  -whorls {NW} \n" \
  "    This mandatory parameter specifies the number of whorls in the stirring.\n" \
  "\n" \
  "  -detailSize {DSIZE} \n" \
  "    Specifies the typical size (in pixels) of the curls created by the rough" \
  " stirring.  The default is " stringify(DEFAULT_DETAIL_SIZE) " pixels.\n" \
  "\n" \
  "  -maxval {MAXVAL} \n" \
  "    Defines the nominal maximum sample value for the output" \
  " image.  Must be an integer between 1" \
  " and " stringify(PNM_MAX_MAXVAL) ".  If not specified," \
  " defaults to " stringify(PNM_MAX_MAXVAL) ".\n" \
  "\n" \
  "  -verbose \n" \
  "    If this flag is present, several diagnostic messages and mask statistics" \
  " are written to {stderr}.  The default is to supress such" \
  " diagnostics.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  pgmkernel(1), pamgauss(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 2012-02-18 by J. Stolfi (IC-UNICAMP).\n" \
  "\n" \
  "HISTORY\n" \
  "  2012-02-18 by J.Stolfi: created.\n" \
  "  2012-02-25 by J.Stolfi: added \"-whorls\" option.\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " ppmmarpap_C_COPY ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define stringify(x) strngf(x)
#define strngf(x) #x

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <values.h>

#include <r2.h> 
#include <r2x2.h> 
#include <r2_extra.h> 
#include <r3.h> 
#include <ix_reduce.h> 
#include <bool.h> 

#include <jspnm.h> 
#include <jsfile.h> 
#include <jsrandom.h> 
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>
#include <uint16_image_read_pnm.h>
#include <float_image.h>
#include <float_image_transform.h>
#include <float_image_paint.h>
#include <float_image_from_uint16_image.h>
#include <float_image_to_uint16_image.h>
#include <frgb_ops.h>
#include <affirm.h> 
#include <argparser.h> 

/* DATA TYPES */

#define INF INFINITY

typedef struct options_t 
  { int NX, NY;          /* Output image width and height. */
    int whorls;              /* Numer of whorls. */
    int detailSize;      /* Size of smallest whorls. */
    char *source;        /* Source image (or NULL to make up stripes). */
    uint16_t maxval; /* Output maxval. */
    bool_t verbose;      /* TRUE prints diagnostics and statistics. */
  } options_t;
  /* Arguments parsed from the command line. */

#define MAX_IMAGE_SIZE (256*1024)
  /* Max image width or height, only for consistency checking. */
  
#define DEFAULT_DETAIL_SIZE 24  

/* INTERNAL PROTOTYPES */

int main(int argc, char* argv[]);
options_t *parse_options(int argc, char **argv);
float_image_t *read_image(char *fname, bool_t verbose);
float_image_t *make_source_image(int NC, int NX, int NY, int ND, bool_t verbose);
float_image_t *make_test_image(int NC, int NX, int NY, int ND, bool_t verbose);
void paint_stripe(float_image_t *img, int ylo, int yhi, int ND, float pix[]);
float_image_t *stir_paint(float_image_t *iimg, int NX, int NY, int NW, int ND, bool_t verbose);
void pick_nice_color(double H, int NC, float pix[]);
void write_image(FILE *wr, float_image_t *oimg, int maxval, bool_t verbose);

void define_stir_parameters(int NX, int NY, int NW, int ND, r2_t wct[], double wrd[], double wan[], bool_t verbose);
  /* Defines the stirring paramters.  The stirring is a series
    of {NW} whirlpools with random centers in the domain {[0_NX]×[0_NY]}, decreasing radii, and random and decreasing angles.
    Whirlpool number {i} has center {wct[i]}, radius {wrd[i]}, and angle {wan[i]}.
    The larger whirlpool has radius comparable to the image size;
    the smaller ones have radius {~ND/2}. */
    
/* IMPLEMENTATIONS */

int main(int argc, char* argv[])
  { /* Command line arguments: */
    options_t *o = parse_options(argc, argv);
    bool_t debug = FALSE;

    int NC;
    float_image_t *iimg = NULL;
    if ((o->source != NULL) && (strlen(o->source) != 0))
      { iimg = read_image(o->source, o->verbose);
        NC = (int)iimg->sz[0];
      }
    else
      { NC = 3;
        if (debug)
          { iimg = make_test_image(NC, o->NX, o->NY, o->detailSize, o->verbose); }
        else
          { iimg = make_source_image(NC, o->NX, o->NY, o->detailSize, o->verbose); }
      }
      
    if (o->verbose)
      { char *fname = jsprintf("/tmp/init.%s", (NC == 3 ? "ppm" : "pgm"));
        FILE *tmp = open_write(fname, TRUE);
        write_image(tmp, iimg, o->maxval, o->verbose);
        fclose(tmp);
        free(fname);
      }
    
    float_image_t *oimg = stir_paint(iimg, o->NX, o->NY, o->whorls, o->detailSize, o->verbose);
    
    write_image(stdout, oimg, o->maxval, o->verbose);
    
    float_image_free(iimg);
    float_image_free(oimg);
    return 0;
  }
    
void write_image(FILE *wr, float_image_t *fimg, int maxval, bool_t verbose)
  {
    int NC = (int)fimg->sz[0];
    bool_t isMask = FALSE;
    uint16_image_t *pimg = float_image_to_uint16_image
      ( fimg, isMask, NC, NULL, NULL, NULL, (uint16_t)maxval, TRUE, verbose );
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pimg, forceplain, verbose);
    fflush(wr);
    uint16_image_free(pimg);
  }

float_image_t *read_image(char *fname, bool_t verbose)
  { 
    FILE *rd = open_read(fname, verbose);
    uint16_image_t *pim = uint16_image_read_pnm_file(rd);
    bool_t isMask = FALSE; /* Assume smooth distr of pixel values. */
    float_image_t *fim = float_image_from_uint16_image(pim, isMask, NULL, NULL, TRUE, verbose);
    uint16_image_free(pim);
    if (rd != stdin) { fclose(rd); }
    return fim;
  } 
  
float_image_t *make_source_image(int NC, int NX, int NY, int ND, bool_t verbose)
  {
    float_image_t *fimg = float_image_new(NC, NX, NY);
    
    /* Generate some horizontal color stripes: */
    float pix[NC];  /* Current stripe color. */
    int ylo = 0;    /* Last Y of current stripe. */
    int NYS_min = (ND+3)/2;
    int NYS_max = (3*ND+5)/2;
    double H = drandom();      /* Stripe hue. */
    while (ylo < NY)
      { pick_nice_color(H, NC, pix);
        int dy =  int32_abrandom(NYS_min, NYS_max); /* Width of next stripe. */
        int yhi = (int)imin(ylo + dy, NY-1);
        if (NY - 1 - yhi < NYS_min) { yhi = NY - 1; }
        paint_stripe(fimg, ylo, yhi, ND, pix);
        ylo = yhi + 1;
        H = H + (0.06*drandom() - 0.03);
        if (drandom() < 0.200) { H = H + 0.5; }
        while (H >= 1.0) { H = H-1; }
        while (H < 0.0) { H = H+1; }
      }
    return fimg;
  }
  
float_image_t *make_test_image(int NC, int NX, int NY, int ND, bool_t verbose)
  {
    float_image_t *fimg = float_image_new(NC, NX, NY);

    /* Fill image with uniform color: */
    float_image_fill(fimg, 0.90f);
    
    /* Draw grid: */
    float pix[NC];  /* Grid line color. */
    int c, x, y;
    for (c = 0; c < NC; c++) { pix[c] = (float)((c + 0.5)/((double)NC)); }
    int NS = 8;
    for (y = 0; y < NY; y += NS)
      { for (x = 0; x < NX; x++)
          { float_image_set_pixel(fimg, x, y, pix); }
      }
    for (y = 0; y < NY; y ++)
      { for (x = 0; x < NX; x += NS)
          { float_image_set_pixel(fimg, x, y, pix); }
      }
    return fimg;
  }

void paint_stripe(float_image_t *img, int ylo, int yhi, int ND, float pix[])
  {
    int NC = (int)img->sz[0];
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];
    demand((ylo >= 0) && (ylo <= yhi) && (yhi < NY), "bad Y range");
    
    int NYS = yhi - ylo + 1;  /* Number of rows in stripe. */
    
    int x, y;
    for (y = ylo; y <= yhi; y++)
      {
        for (x = 0; x < NX; x++)  { float_image_set_pixel(img, x, y, pix); }
      }
      
    double rmax = 0.07*((double)ND); /* Max radius of black spots. */
    double dmin = 0.02*((double)ND); /* Min distance from spot to band edge. */
    double hwd = 0.0; /* Dot outline width (none). */
    float vdraw = NAN; /* Do not draw outline. */
    float vfill = 0.0; /* Fill with black. */
    
    if (NYS >= 2*(rmax+dmin))
      { 
        /* Add some random black spots: */
        int NG = (int)ceil(0.125*((double)NX*NYS)/((double)ND*ND)); /* Num of spot groups. */
        bool_t round = TRUE;
        bool_t diagonal = FALSE;
        int sub = 4;
        int kg;
        for (kg = 0; kg < NG; kg++)
          { double xctr = drandom()*NX;
            double yctr = dmin + rmax + (NYS - 2*(dmin + rmax))*drandom(); /* Rel to {ylo}. */
            /* Num of spots in group: */
            int NS = abs(int32_abrandom(-3,+3) + int32_abrandom(-3,+3) + int32_abrandom(-3,+3)) + 1;  
            int ks;
            for (ks = 0; ks < NS; ks++)
              { double rad = (0.60 + 0.40*drandom())*rmax;
                int c;
                for (c = 0; c < NC; c++)
                  { (void)float_image_paint_dot
                      ( img, c, xctr, ylo + yctr, rad, hwd, 
                        round, diagonal, vfill,vdraw, sub
                      );
                  }
                xctr = xctr + (0.45 + 0.15*drandom())*((double)ND);
                yctr = yctr + (0.06*drandom() - 0.03)*((double)ND);
                yctr = fmax(yctr, dmin + rmax);
                yctr = fmin(yctr, NYS - dmin - rmax);
              }
          }
      }
  }

float_image_t *stir_paint(float_image_t *iimg, int NX, int NY, int NW, int ND, bool_t verbose)
  {
    bool_t debug = FALSE;
    
    int NC =  (int)iimg->sz[0];
    int NXI = (int)iimg->sz[1];
    int NYI = (int)iimg->sz[2];
    float_image_t *oimg = float_image_new(NC, NX, NY);
    ix_reduce_mode_t red = ix_reduce_mode_MIRROR;
    bool_t avg = TRUE;
    int order = 0; /* Interpolation order */

    /* Stirring parameters: */
    r2_t wct[NW];    /* Whirlpool centers. */
    double wrd[NW];  /* Whirlpool centers. */
    double wan[NW];  /* Whirlpool angles. */
    define_stir_parameters(NXI, NYI, NW, ND, wct, wrd, wan, verbose);
    
    auto void stir_map (r2_t *p, r2x2_t *J);
      /* Stores in {q} the point of {[0_NX]×[0_NY]}
        that ends up at {p} after the stirring .*/
    
    float_image_transform_all(iimg, red, &stir_map, 0.0, avg, order, NULL, oimg); 
    
    if (debug)
      { bool_t round = TRUE;
        bool_t diagonal = FALSE;
        int sub = 4;
        double drad = 2.0; /* Dot radius. */
        double dhwd = 0.5; /* Pen half-width. */
        int i;
        for (i = 0; i < NW; i++) 
          { int c;
            r2_t ctr = wct[i];
            double rad = wrd[i];
            for (c = 0; c < NC; c++)
              { (void)float_image_paint_dot(oimg, c, ctr.c[0], ctr.c[1], drad,   0, round, diagonal, 1.0, NAN, sub); 
                (void)float_image_paint_dot(oimg, c, ctr.c[0], ctr.c[1], rad, dhwd, round, diagonal, NAN, 1.0, sub);
              }
          }
      }
      
    return oimg;
    
    void stir_map (r2_t *p, r2x2_t *J)
      {
        /* Map domain of {oimg} to the domain of {iimg}: */
        double sx = ((double)NXI)/((double)NX);
        double sy = ((double)NYI)/((double)NY);
        p->c[0] = sx*p->c[0];
        p->c[1] = sy*p->c[1];
        if (J != NULL) 
          { r2x2_t PJ;
            r2x2_ident(&PJ);
            PJ.c[0][0] = sx;
            PJ.c[1][1] = sy;
            r2x2_mul(J, &PJ, J);
          }

        int i;
        for (i = 0; i < NW; i++) { r2_map_twirl(p, &(wct[i]), wrd[i], wan[i], J); }
      }
  }

void define_stir_parameters(int NX, int NY, int NW, int ND, r2_t wct[], double wrd[], double wan[], bool_t verbose)
  { 
    double radMax = 0.20*fmin(NX, NY);        /* Max whorl radius. */
    double radIni = 0.20*hypot(NX, NY);       /* Ideal radius of first whorl. */
    double radFin = fmin(2*ND,0.50*radIni);   /* Ideal radius of last whorl. */
    double scaFin = radFin/radIni;            /* Ideal scale from last to first. */
    assert(scaFin < 1.0);
    /* Coeffs of scale formula: */
    double beta = NW;
    double alfa = (scaFin*sqrt(beta+NW-1) - sqrt(beta))/(1 - scaFin);
    if (verbose) 
      { fprintf(stderr, "  radMax = %14.8f  radIni = %14.8f  radFin = %14.8f\n", radMax, radIni, radFin);
        fprintf(stderr, "  scaFin = %14.8f  alfa =   %+14.8f  beta =   %+14.8f\n", scaFin, alfa, beta);
      }
    int i, k;
    for (i = 0; i < NW; i++)
      { /* Decide the relative scale of this whorl: */
        double scale = (alfa + sqrt(beta))/(alfa + sqrt(i + beta));
        /* Pick a turn angle {ang} in radians: */
        double ang = 0;
        for (k = 0; k < 4; k++)
          { double ank = sqrt(scale)*3*M_PI*(drandom()-drandom());
            if (fabs(ank) > fabs(ang)) { ang = ank; }
          }
        /* Pick a nominal radius {rad} in image space: */
        double rad = fmin(radMax, scale*radIni);
        r2_t ctr;
        ctr = (r2_t){{ NX*drandom(),NY*drandom() }};
        if (verbose)
          { r2_gen_print(stderr, &ctr, "%14.8f", "ctr = ( ", " ", " )");
            fprintf
              ( stderr, "  scale = %14.8f  rad = %14.8f  ang = %+12.8f (%+8.2f°)\n", 
                scale, rad, ang, ang*180/M_PI
              );
          }
        
        /* Set param vector: */
        wct[i] = ctr;
        wrd[i] = rad;
        wan[i] = ang;
      }
  }

void pick_nice_color(double H, int NC, float pix[])
  { 
    if (NC == 1)
      { double V = drandom();
        pix[0] = (float)V;
      }
    else
      { double V = (drandom() < 0.300 ? 0.3 : 0.8 + 0.2*drandom());
        double S = V;
        frgb_t color = (frgb_t){{ (float)H, (float)S, (float)V }};
        frgb_print(stderr, " HSV = ( ", &color, 3, "%5.3f", " ) ");
        frgb_from_HSV_CG(&color);
        frgb_print(stderr, " RGB = ( ", &color, 3, "%5.3f", " )\n");
        int c;
        for (c = 0; c < NC; c++) { pix[c] = color.c[c % 3];}
      }
  }

options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 
    
    /* Get keyword arguments: */
    argparser_get_keyword(pp, "-size");
    o->NX = (int)argparser_get_next_int(pp, 1, MAX_IMAGE_SIZE); 
    o->NY = (int)argparser_get_next_int(pp, 1, MAX_IMAGE_SIZE); 
    
    if (argparser_keyword_present(pp, "-source"))
      { o->source = argparser_get_next_non_keyword(pp); }
    else
      { o->source = NULL; }
      
    argparser_get_keyword(pp, "-whorls");
    o->whorls = (int)argparser_get_next_int(pp, 0, 10000);
        
    if (argparser_keyword_present(pp, "-detailSize"))
      { o->detailSize = (int)argparser_get_next_int(pp, 0, MAX_IMAGE_SIZE); }
    else
      { /* Use input maxval: */
        o->detailSize = DEFAULT_DETAIL_SIZE;
      }
    
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (uint16_t)argparser_get_next_int(pp, 0, PNM_FILE_MAX_MAXVAL); }
    else
      { /* Use input maxval: */
        o->maxval = PNM_FILE_MAX_MAXVAL;
      }
    
    /* Diagnostics and statistics: */
    o->verbose = argparser_keyword_present(pp, "-verbose");
    
    /* Skip to positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }
