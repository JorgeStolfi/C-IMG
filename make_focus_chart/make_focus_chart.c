#define PROG_NAME "make_focus_chart"
#define PROG_DESC "creates an EPS file with a simple focus calibration chart"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-21 13:59:58 by stolfi */

#define make_focus_chart_C_COPYRIGHT \
  "Copyright © 2017  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -size {NX} {NY} \\\n" \
  "    -spacing {SPACING} \\\n" \
  "    [ -radius {RADIUS} ] \\\n" \
  "    [ -lineWidth {LINEWIDTH} ] \\\n" \
  "    [ -nCopies {NCPX} {NCPY} ] \\\n" \
  "    [ -dpi {DPI} ] \\\n" \
  "    [ -geomFile {GFNAME} ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Writes to standard ouput a EPS file with the drawing" \
  " of a /focus chart/ a small patterned chart" \
  " that can be used to check focus and focal depth of a low-power" \
  " microscope.\n" \
  "\n" \
  "  The chart consists of a square grid of black crosses and open" \
  " circles of the given radius and line width" \
  " (in mm) against a white background.  The crosses" \
  " are rotated by random angles.\n" \
  "\n" \
  "  Optionally writes also a text description of the chart.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -size {NX} {NY}\n" \
  "    Specifies the number of grid marks along each axis.\n" \
  "\n" \
  "  -spacing {SPACING}\n" \
  "    This mandatory parameter specifies the horizontal and vertical" \
  " spacing between the centers of the crosses or circles.  It will be" \
  " rounded to an even multiple of the dot size.\n" \
  "\n" \
  "  -radius {RADIUS}\n" \
  "    This optional parameter specifies the radius of each cross or" \
  " circle in millimeters. If not specified, it defaults to {SPACING/3}.  It" \
  " will berounded to an even multiple of the dot size.\n" \
  "\n" \
  "  -lineWidth {LINEWIDTH}\n" \
  "    This optional parameter specifies the width (in millimeters) of the lines used" \
  " to draw the crosses. If not specified, it defaults to {SPACING/6}. It" \
  " will be rounded to an even multiple of the dot size\n" \
  "\n" \
  "  -nCopies {NCPX} {NCPY}\n" \
  "    This optional parameter specifies the number of copies of" \
  " the chart to draw on the page. They will be arranged" \
  " into {NCPX} columns and {NCPY} rows.\n" \
  "\n" \
  "  -dpi {DPI}\n" \
  "    This optional parameter specifies the printer resolution to be" \
  " assumed, in dots per inch. All dimensions will be rounded" \
  " to a multiple of this parameter.  If not specified, defaults to 600.\n" \
  "\n" \
  "  -geomFile {GFNAME}\n" \
  "    This optional parameter specifies the name of the file which will" \
  " contain the coordinates of centers and rotation angles of the crosses.  If omitted," \
  " this information is not written.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgmtopsx(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2017-06-07 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2017-06-07 Created from {make_light_chart.c}. J.Stolfi.\n" \
  "  2018-07-02 Separated the floor grid to {make_mark_grid.c}. J.Stolfi.\n" \
  "  2018-07-02 Added the \"-nCopies\" option. J.Stolfi.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " make_focus_chart_C_COPYRIGHT "\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <i2.h>
#include <jsfile.h>
#include <epswr.h>
#include <frgb.h>
#include <argparser.h>

#include <mkgr_mark.h>
#include <mkgr_mark_grid.h>
#include <mkgr_mark_grid_draw_eps.h>

typedef struct options_t
  { /* Parameters of focus chart: */
    i2_t size;        /* Num columns and num rows of marks. */
    double spacing;   /* Spacing (mm) of marks along each axis. */
    double radius;    /* Radius of crosses (mm). */
    double lineWidth; /* Line width (mm) for crosses and open circles. */
    char *geomFile;   /* Output file name with marks data. */
    /* General parameters: */
    i2_t nCopies;     /* Number of copies to print on the page. */
    double dpi;       /* Printer resolution to assume. */
  } options_t;
  
/* PROTOTYPES */

int main (int argc, char **argv);

mkgr_mark_grid_t *mfc_create_focus_chart_grid(int NX, int NY, double spacing, double rcross, double lwd, double dpi);
  /* Creates the mark list for the focus chart proper, consisting of
    {NX} columns and {NY} rows of marks, with the given {spacing} (mm)
    beteen centers. The marks will be crosses of radius {rcross} (mm),
    in various orientations, and circles of radius proportional to
    {rcross}.  Crosses and open circles will be drawn with pen of
    width {lws} (mm). Assumes that the printer has {dpi} dots per inch. */
    
void mfc_draw_figure(options_t *o, mkgr_mark_grid_t *gr);
  /* Draws the figure, consisting of multiple copies of the marks in the 
    grid {gr}. */
  
void mfc_draw_focus_chart(epswr_figure_t *eps, mkgr_mark_grid_t *gr, r2_t pos, r2_t tSize, double dpi);
  /* Draws a focus chart.  
    
    The chart consists of four sections strung side by side. Each
    section consists of a rectangle of size {tSize}.
    Each of the two middle sections contains a copy of the grid {gr},
    centered in it. */

void mfc_get_size_focus_chart(mkgr_mark_grid_t *gr, double mrg, r2_t *sMin, r2_t *sMax, r2_t *tSize);
  /* Returns in {*sMin,*sMax} the bounding box of the focus chart,
    relative to the low corner of the first section, incluing the marks,
    frames, and gluing tabs, with a margin {4*mrg} all around. Also
    returns in {*tSize} the width and height of each section, which are
    the width and height of the mark grid plus an inner margin {mrg}
    around it. */

options_t *mfc_parse_options(int argc, char **argv);
  /* Obtains arguments from the command line. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  { 
    /* Command-line parameters: */
    options_t *o = mfc_parse_options(argc, argv);
    
    /* Mark grid parameters (unrounded): */
    double dpi = o->dpi;
    double spc = o->spacing;
    double radius = o->radius;
    double lwd = o->lineWidth;
    int NX = o->size.c[0];
    int NY = o->size.c[1];
    mkgr_mark_grid_t *gr = mfc_create_focus_chart_grid(NX, NY, spc, radius, lwd, dpi);
    
    mfc_draw_figure(o, gr);
    if(o->geomFile != NULL)
      { /* Write chart geometry file: */
        FILE *geomFile = open_write(o->geomFile, TRUE);
        mkgr_mark_grid_write(geomFile, gr);
        fclose(geomFile);
      }
    mkgr_mark_grid_free(gr);  
    return 0;
  }

mkgr_mark_grid_t *mfc_create_focus_chart_grid(int NX, int NY, double spc, double radius, double lwd, double dpi)  
  { fprintf(stderr, "nominal: spacing = %.5f  radius = %.5f  lineWidth = %.5f\n", spc, radius, lwd); 
    spc = epswr_round_dim_mm(spc, dpi);
    lwd = epswr_round_dim_mm(lwd, dpi);
    double rcross = epswr_round_dim_mm(radius, dpi);
    double rcircle = epswr_round_dim_mm(radius, dpi); /* Circle radius in mm. */
    double rdot = epswr_round_dim_mm(lwd/2, dpi);
    fprintf(stderr, "actual: spacing = %.5f", spc);
    fprintf(stderr, "  rcross = %.5f  rcircle = %.5f  rdot = %.5f", rcross, rcircle, rdot);
    fprintf(stderr, "  lineWidth = %.5f\n", lwd); 

    auto mkgr_mark_t defmark_main(int ix, int iy);
      /* Crosses alternating upright and diagonal, with circles in quincux of step 4. */
     
    auto mkgr_mark_t defmark_second(int ix, int iy);
      /* Small dots in the middle of cells of main grid. */
    
    mkgr_mark_grid_t *gr = mkgr_mark_grid_new();
    
    i2_t iMin = (i2_t){{ 0, 0 }};
    i2_t iMax = (i2_t){{ NX-1, NY-1 }};
    mkgr_mark_grid_append_marks(gr, iMin, iMax, defmark_main);
    
    i2_t iMin2 = (i2_t){{ 0, 0 }};
    i2_t iMax2 = (i2_t){{ NX-2, NY-2 }};
    mkgr_mark_grid_append_marks(gr, iMin2, iMax2, defmark_second);
    return gr;
    
    mkgr_mark_t defmark_main(int ix, int iy)
      { mkgr_mark_t mk;
        double xctr = ix*spc; /* Assumes spc is rounded to an even multiple of dot size. */
        double yctr = iy*spc; /* Assumes spc is rounded to an even multiple of dot size. */
        mk.ctr = (r2_t){{ xctr, yctr }};
        bool_t cross = ((iy % 2) != 0) || (((ix+iy) % 4) != 0);
        mk.cross = cross;
        mk.rad = (cross ? rcross : rcircle);
        mk.color = (frgb_t){{ 0.000, 0.000, 0.000 }};
        mk.lwd = lwd;
        bool_t diag = ((ix+iy) % 2 == 1);
        double ang = (diag & cross ? 0.125 : 0.000); 
        mk.ang = ang - floor(ang);
        return mk;
      }
    
    mkgr_mark_t defmark_second(int ix, int iy)
      { mkgr_mark_t mk;
        double xctr = (ix+0.5)*spc; /* Assumes spc is rounded to an even multiple of dot size. */
        double yctr = (iy+0.5)*spc; /* Assumes spc is rounded to an even multiple of dot size. */
        mk.ctr = (r2_t){{ xctr, yctr }};
        mk.cross = FALSE;
        mk.rad = rdot;
        mk.color = (frgb_t){{ 0.000, 0.000, 0.000 }};
        mk.lwd = 0.0;
        mk.ang = 0.0;
        return mk;
      }
  }

void mfc_draw_figure(options_t *o, mkgr_mark_grid_t *gr)
  { 
    /* The figure has {o->nCopies} copies of the focus charts. */
    /* Each focus chart has two copies of the mark grid, and two gluing tabs, along X axis. */
    
    bool_t debug = FALSE;
   
    /* Get bounding box of the focus chart: */
    double mrg = epswr_round_dim_mm(0.16, o->dpi); /* Extra margin around each mark grid and whole chart. */
    r2_t pMin, pMax; /* Low and high corners of focus chart, rel to its origin, including margins (mm). */
    r2_t tSize; /* Width and height of each section of the focus chart. */ 
    mfc_get_size_focus_chart(gr, mrg, &pMin, &pMax, &tSize);
    r2_t pSize; r2_sub(&pMax, &pMin, &pSize); /* Width and height of the focus chart, with all margins. */
    fprintf(stderr, "chart size ( %9.6f %9.6f )", pSize.c[0], pSize.c[1]);
    fprintf(stderr, "  section size ( %9.6f %9.6f )\n", tSize.c[0], tSize.c[1]);
    
    /* Position the focus charts: */
    int NCX = o->nCopies.c[0]; /* Number of columns of copies of the chart. */
    int NCY = o->nCopies.c[1]; /* Number of rows of copies of the chart. */
    int NC = NCX*NCY; /* Total number of copies of the chart. */
    r2_t org[NC]; /* Origin of each focus chart, rel to {fMin} (mm). */
    int kc = 0;
    for (int kx = 0; kx < NCX; kx++)
      { for (int ky = 0; ky < NCY; ky++)
         { /* Compute the lower corner of this copy: */
           r2_t pLo = (r2_t){{ kx*pSize.c[0], ky*pSize.c[1] }};
           /* Compute the origin of this copy: */
           r2_t orgk; r2_sub(&pLo, &pMin, &(orgk));
           /* For good measure, round the coordinates: */
           for (int j = 0; j < 2; j++) {  orgk.c[j] = epswr_round_dim_mm(orgk.c[j], o->dpi); }
           if (debug) { fprintf(stderr, "placing chart at ( %9.6f %9.6f )\n", orgk.c[0], orgk.c[1]); }
           org[kc] = orgk;
           kc++;
         }
      }
    assert(kc == NC);
    
    /* Figure bounding box (mm), relative to origin of figure: */
    r2_t fMin = (r2_t){{ 0.0, 0.0 }}; /* Low corner. */
    r2_t fMax = (r2_t){{ NCX*pSize.c[0], NCY*pSize.c[1] }};; /* High corner. */
    
    /* Total size of plot area: */
    double us_szx_mm = fMax.c[0] - fMin.c[0];
    double us_szy_mm = fMax.c[1] - fMin.c[1];
    
    /* Total size of figure (including margins outside of plot area, if {eps} is true): */
    double fig_mrg_mm = epswr_round_dim_mm(2.0, o->dpi), fig_mrg_pt = fig_mrg_mm * epswr_pt_per_mm;
    double fig_szx_mm = us_szx_mm + 2*fig_mrg_mm, fig_szx_pt = fig_szx_mm * epswr_pt_per_mm;
    double fig_szy_mm = us_szy_mm + 2*fig_mrg_mm, fig_szy_pt = fig_szy_mm * epswr_pt_per_mm;
    
    /* Open the figure stream/document: */
    bool_t verbose = TRUE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( "out", NULL, "chart", -1, NULL,
        fig_szx_pt, fig_szy_pt, 
        fig_mrg_pt, fig_mrg_pt, fig_mrg_pt, fig_mrg_pt,
        verbose
      );
    epswr_set_label_font(eps,"Courier",10.0);
    epswr_set_pen(eps, 1.000, 0.000, 0.000,  0.40,  0.0, 0.0);
    
    /* Set up the plot window for the whole plot: */
    double hMin_pt = fig_mrg_mm*epswr_pt_per_mm, hMax_pt = hMin_pt + us_szx_mm*epswr_pt_per_mm;
    double vMin_pt = fig_mrg_mm*epswr_pt_per_mm, vMax_pt = vMin_pt + us_szy_mm*epswr_pt_per_mm;
    double xMin_mm = fMin.c[0], xMax_mm = fMax.c[0];
    double yMin_mm = fMin.c[1], yMax_mm = fMax.c[1];
    epswr_set_window(eps, hMin_pt, hMax_pt, vMin_pt, vMax_pt, FALSE,  xMin_mm, xMax_mm, yMin_mm, yMax_mm);
    
    /* Draw the focus charts: */
    kc = 0;
    for (int kx = 0; kx < NCX; kx++)
      { for (int ky = 0; ky < NCY; ky++)
         { mfc_draw_focus_chart(eps, gr, org[kc], tSize, o->dpi);
           kc++;
        }
      }
    assert(kc == NC);
    epswr_end_figure(eps);
  }

void mfc_draw_focus_chart(epswr_figure_t *eps, mkgr_mark_grid_t *gr, r2_t org, r2_t tSize, double dpi)
  { 
    double lineWidth_frame = 0.1;
    
    /* Daw the four sections: */
    double xlo = org.c[0]; /* Low X of tab incl. margin */
    double ylo = org.c[1]; /* Low Y of tab incl. margin */
    for (int kt = 0; kt < 4; kt++)
      { /* Draw frame around section {kt}: */
        double xhi = xlo + tSize.c[0];
        double yhi = ylo + tSize.c[1];
        epswr_set_pen(eps, 0.000, 0.000, 0.000,  lineWidth_frame,  0.0, 0.0);
        epswr_set_fill_color(eps, NAN, NAN, NAN);
        epswr_rectangle(eps, xlo, xhi, ylo, yhi, FALSE, TRUE);
        if ((kt == 1) || (kt == 2))
          { /* Compute the position {(gx,xy)} of the grid's origin: */
            double gx = 0.5*((xlo + xhi) - (gr->pMax.c[0] + gr->pMin.c[0]));
            double gy = 0.5*((ylo + yhi) - (gr->pMax.c[1] + gr->pMin.c[1]));
            r2_t gorg = (r2_t){{ gx, gy }};
            /* Draw the mark grid inside section: */
            mkgr_mark_grid_draw_eps(eps, gr, &gorg, dpi);
          }
        /* Move to next tab: */
        xlo += tSize.c[0];
      }
  }
         
void mfc_get_size_focus_chart(mkgr_mark_grid_t *gr, double mrg, r2_t *sMin, r2_t *sMax, r2_t *tSize)
  { 
    /* Height and width of each section: */
    double sx = gr->pMax.c[0] - gr->pMin.c[0] + 2*mrg;
    double sy = gr->pMax.c[1] - gr->pMin.c[1] + 2*mrg;
    (*tSize) = (r2_t){{ sx, sy }};
    
    /* Origin is at low corner of first section: */
    (*sMin) = (r2_t){{ -4*mrg, -4*mrg }};
    (*sMax) = (r2_t){{ 4*sx + 4*mrg, sy + 4*mrg }};
  }

options_t *mfc_parse_options(int argc, char **argv)
  {
    /* Initialize the argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    argparser_get_keyword(pp, "-size");
    o->size.c[0] = (int32_t)argparser_get_next_int(pp, 0, 1000);
    o->size.c[1] = (int32_t)argparser_get_next_int(pp, 0, 1000);

    argparser_get_keyword(pp, "-spacing");
    o->spacing = argparser_get_next_double(pp, 0, 1000.0);

    if (argparser_keyword_present(pp, "-dpi"))
      { o->dpi = argparser_get_next_double(pp, 36, 1000.0); }
    else
      { o->dpi = 600.0; }

    if (argparser_keyword_present(pp, "-radius"))
      { o->radius = argparser_get_next_double(pp, 0, 1000.0);  }
    else
      { o->radius = o->spacing/3; }

    if (argparser_keyword_present(pp, "-lineWidth"))
      { o->lineWidth = argparser_get_next_double(pp, 0, 1000.0);  }
    else
      { o->lineWidth = o->spacing/6; }

    if (argparser_keyword_present(pp, "-nCopies"))
      { o->nCopies.c[0] = (int32_t)argparser_get_next_int(pp, 1, 5);
        o->nCopies.c[1] = (int32_t)argparser_get_next_int(pp, 1, 50);
      }
    else
      { o->nCopies = (i2_t){{ 2, 3 }}; }

    if (argparser_keyword_present(pp, "-geomFile"))
      { o->geomFile = argparser_get_next(pp); }
    else
      { o->geomFile = NULL; }

    /* Get positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for leftover args: */
    argparser_finish(pp);

    return o;
  }

