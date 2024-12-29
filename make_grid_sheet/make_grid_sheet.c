#define PROG_NAME "make_grid_sheet"
#define PROG_DESC "creates a Postscript and/or PNG file with a grid of dots and crosses"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-21 13:59:56 by stolfi */

#define make_grid_sheet_C_COPYRIGHT \
  "Copyright © 2018  by the State University of Campinas (UNICAMP)"
  
/* !!! TO DO: Option to invert colors of marks and background. !!! */
/* !!! TO DO: Make disk optional. !!! */
/* !!! TO DO: Option to rotate marks by 'random' angles. !!! */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -gridSize {NX} {NY} \\\n" \
  "    [ -markSpacing {MARK_SPACING} ] \\\n" \
  "    [ -markRadius {MARK_RADIUS} ] \\\n" \
  "    [ -lineWidth {LINE_WIDTH} ] \\\n" \
  "    [ -printerDPI {PRINTER_DPI} ] \\\n" \
  "    [ -imagePXPMM {IMAGE_PXPMM} ] \\\n" \
  "    -outPrefix {OUT_PREFFIX}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Writes to disk a Postscript file \"{OUTPREFIX}.eps\" and/or a PNG" \
  " image \"{OUTPREFIX}.png\" with a grid of marks (dots, circles, or crosses)" \
  " that can be used as a background when photographing small objects.\n" \
  "\n" \
  "  The level 1 marks will be open circles with center dots, at the" \
  " corners of all 5x5 blocks of cells.  The level 2 marks are" \
  " upright crosses at all cell corners, excluding the" \
  " level 1 ones.  The level 3 marks will be small dots" \
  " at the center of the cells, except when too close" \
  " to a level 1 mark.\n" \
  "\n" \
  "  The Postscript version of the grid consists of a square grid of or white crosses, dots, and open" \
  " circles, against a black background disk.  Mark centers will be rounded to an" \
  " integer number of printer dots, but mark spacings may vary by 1 dot.  Marks that lie" \
  " outside the disk will be invisible.\n" \
  "\n" \
  "  The PNG grid has the marks against a full black background.  The" \
  " mark spacing will be rounded to an even integer number of pixels.  All marks" \
  " will be centered on pixel corners.a" \
  " The image will span only" \
  " the mark centers, exactly. Therefore, the outermost marks will be clipped in half.\n" \
  "\n" \
  "  Also writes text files \"{OUTPREFIX}_eps.txt\" and/or \"{OUTPREFIX}_png.txt\"containing" \
  " the parameters of all marks (center coordinates, radius, mark" \
  " type, rotation angle, etc.) in the corresponding graphics file.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -gridSize {NX} {NY}\n" \
  "    Specifies the number of grid cells along each axis. Must be even.\n" \
  "\n" \
  "  -markRadius {MARK_RADIUS}\n" \
  "    Specifies the radius of each cross or circle in millimeters.  Currently" \
  " applies only to the EPS version.\n" \
  "\n" \
  "  -lineWidth {LINE_WIDTH}\n" \
  "    Specifies the width (in millimeters) of the lines used" \
  " to draw the crosses and circle outlines.\n" \
  "\n" \
  "  -markSpacing {MARK_SPACING}\n" \
  "    This optional parameter specifies the horizontal and vertical" \
  " spacing between the centers of the main grid marks, in mm. If not" \
  " specified, it defaults to 1 mm.\n" \
  "\n" \
  "  -printerDPI {PRINTER_DPI}\n" \
  "    This optional parameter specifies the printer resolution to be" \
  " assumed when generating the Postscript file, in dots" \
  " per inch.  Mark dimensions and the line width will be rounded" \
  " to a multiple of this parameter.  If not specified or zero, the Postscript" \
  " file will not be written.\n" \
  "\n" \
  "  -imagePXPMM {IMAGE_PXPMM}\n" \
  "    This optional parameter specifies the resolution to be assumed for the PNG image," \
  " in pixels per millimeter.  The actual resolution will be rounded up so that" \
  " the mark spacing is an integer number of pixels. If not given, or is zero, the" \
  " PNG image will not be written.\n" \
  "\n" \
  "  -outPrefix {OUTPREFIX}\n" \
  "    This mandatory parameter specifies the prefix for output file names.\n" \
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
  "  2018-07-02 Created from {make_focus_chart.c}. J.Stolfi.\n" \
  "  2018-07-02 Added option to write an image. J.Stolfi.\n" \
  "  2020-11-29 Revamped and switched from {pswr.h} to {epswr.h}. J.Stolfi.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " make_grid_sheet_C_COPYRIGHT "\n" \
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
#include <jsstring.h>
#include <epswr.h>
#include <frgb.h>
#include <argparser.h>

#include <float_image.h>
#include <float_image_write_gen.h>

#include <mkgr_mark.h>
#include <mkgr_mark_grid.h>
#include <mkgr_mark_grid_parms.h>
#include <mkgr_mark_grid_draw_eps.h>
#include <mkgr_mark_grid_paint_image.h>

typedef struct options_t
  { /* Parameters of focus chart: */
    i2_t gridSize;        /* Num columns and num rows of cells (must be even). */
    double markSpacing;   /* Spacing (mm) of marks along each axis. */
    double markRadius;    /* Radius of crosses (mm). */
    double lineWidth;     /* Line width (mm) for crosses and open circles. */
    char *outPrefix;      /* Prefix for output file names. */
    /* General parameters: */
    double printerDPI;    /* Printer resolution in dots per inch. */
    double imagePXPMM;     /* Image resolution in pixels per millimeter. */
  } options_t;
  
/* PROTOTYPES */

int32_t main (int32_t argc, char **argv);

mkgr_mark_grid_t *mgs_create_mark_grid
  ( int32_t HGX, 
    int32_t HGY, 
    mkgr_mark_grid_parms_t *parms
  );
  /* Creates the mark list for a grid with {HGX} by {HGY} cells on each
    side of each axis. There will be a total of {2*HGX+1} by {2*HGY+1}
    white marks, with specified {spacing} between mark centers, mark radii
    {parms.rcross,parms.rcircle,parms.rdot}, and line width {parms.lwd}. 
    
    The main marks will be upright crosses, except that circles will be 
    drawn instead every 5x5 marks, including at the very center. 
    There may be secondary marks besides those main marks.  The 
    marks will be in color {R,G,B}.
    
    If {parms.rbgcirc} is positive, a black background circle with that radius
    will be drawn before all marks. */

void mgs_draw_grid_eps(options_t *o);
  /* Draws the mark grid as an Encapsulated Postscript file
    "{o.outPrefix}.eps", if {o->imageDPI} is not 0.  
    The grid will have white marks and a black background disk of 
    contrasting color. 
    
    Mark sizes and line widths are rounded to integer number of
    printer dots.  The mark spacing however is not rounded.
    
    Also writes a text file "{o.outPrefix}-eps.txt" with the parameters
    of all marks. */
  
void mgs_draw_grid_png(options_t *o);
  /* Draws the mark grid as a PNG image, if {o->imagePXPMM} is not zero.
    Uses white images on a black background.
    
    Rounds {o->imagePXPMM} so that the mark spacing is an integer number
    of pixels. The image will span the mark centers; therefore the
    outermost marks will be clipped. 
    
    Also writes the text file "{o.outPrefix}-png.txt" with the parameters
    of all marks. */
    
void mgs_get_grid_bounding_box(mkgr_mark_grid_t *gr, double mrg, r2_t *sMin, r2_t *sMax);
  /* Returns the bounding box for the grid described by {gr}, including
    all marks, assuming both are centered at the origin; plus a margin
    of {mrg} mm all around.  The box is defined by the low ahd high corners {sMin,sMax}. */
  
void mgs_write_png_image(char *prefix, float_image_t *img);
  /* Writes the image {img} as a PNG file with name "{prefix}.png". */

void mgs_write_mark_data(char *prefix, char *tag, mkgr_mark_grid_t *gr);
  /* Writes the coordinates of the marks in {gr} to a text
    file with name "{prefix}_{tag}.txt". */

options_t *mgs_parse_options(int32_t argc, char **argv);
  /* Obtains arguments from the command line. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { 
    /* Command-line parameters: */
    options_t *o = mgs_parse_options(argc, argv);
    demand(o->gridSize.c[0] % 2 == 0, "grid cols must be even");
    demand(o->gridSize.c[1] % 2 == 0, "grid rows must be even");
    
    /* Spacing of the main marks: */
    
    mgs_draw_grid_eps(o);
    mgs_draw_grid_png(o);
      
    return 0;
  }

void mgs_draw_grid_eps(options_t *o)
  { 
    /* Nominal resolution of medium: */
    double dpi = o->printerDPI; 
    if (dpi <= 0) { return; }

    fprintf(stderr, "--- writing the EPS file ---\n");

    /* Number of grid cells on each side of the axes: */
    int32_t HGX = o->gridSize.c[0]/2;
    int32_t HGY = o->gridSize.c[1]/2;
    
    /* Mark spacing is as  specified, in mm, not rounded: */
    double mkrad = o->markRadius;    /* Requested mark radius (mm), unrounded. */

    /* Define the mark parameters, rounded to integer multiples of printer dots: */
    int32_t HGmin = (int32_t)imin(HGX, HGY);
    double rbgcirc = (HGmin - 0.5)*o->markSpacing;
    mkgr_mark_grid_parms_t parms;
    parms.rcross = 0.75*mkrad;
    parms.rcircle = mkrad;
    parms.rdot = 0.375*mkrad;
    parms.lwd = o->lineWidth;
    parms.rbgcirc = rbgcirc;
    parms.spacing = o->markSpacing;
    mkgr_mark_grid_parms_round(&parms, dpi, 25.4);
    
    /* Create the mark grid {gr}: */
    mkgr_mark_grid_t *gr = mgs_create_mark_grid(HGX, HGY, &parms);
    
    /* Get bounding box of the grid: */
    double mrg_mm = epswr_round_dim_mm(1.0, dpi); /* Margin inside and outside plot area: */
    r2_t pMin, pMax; /* Low and high corners of grid including margin (mm). */
    mgs_get_grid_bounding_box(gr, mrg_mm, &pMin, &pMax);

    /* Size of plot area (pt): */
    r2_t pSize; r2_sub(&pMax, &pMin, &pSize); /* Width and height of the grid, with bg disk and margins. */
    double hPlotSize = pSize.c[0] * epswr_pt_per_mm;
    double vPlotSize = pSize.c[1] * epswr_pt_per_mm;

    /* Open the figure stream/document: */
    char *fname = txtcat(o->outPrefix, ".eps");
    FILE *wr = open_write(fname, TRUE);
    bool_t eps_verb = TRUE;
    double mrg_pt = mrg_mm * epswr_pt_per_mm;
    epswr_figure_t *eps = epswr_new_figure(wr, hPlotSize, vPlotSize, mrg_pt, mrg_pt, mrg_pt, mrg_pt, eps_verb);
    
    epswr_set_label_font(eps,"Courier",10.0);
    
    epswr_set_client_window(eps, pMin.c[0], pMax.c[0], pMin.c[1], pMax.c[1]);
    
    /* Draw the marks: */
    mkgr_mark_grid_draw_eps(eps, gr, NULL, dpi);

    epswr_end_figure(eps);
    
    mgs_write_mark_data(o->outPrefix, "eps", gr);
    
    fprintf(stderr, "--- EPS file written ---\n");

    mkgr_mark_grid_free(gr);
  }

void mgs_draw_grid_png(options_t *o)
  { 
    double pxpmm = o->imagePXPMM; /* Approximate pixels per mm. */
    if (pxpmm <= 0) { return; }
    
    fprintf(stderr, "--- writing the PNG image file ---\n");

    /* Number of grid cells on each side of the axes: */
    int32_t HGX = o->gridSize.c[0]/2;
    int32_t HGY = o->gridSize.c[1]/2;
    
    /* Define the mark parameters, rounded to integer multiples of printer dots: */
    mkgr_mark_grid_parms_t parms;
    parms.rcross = 0.3*o->markSpacing;
    parms.rcircle = 0.4*o->markSpacing;
    parms.rdot = 0.15*o->markSpacing;
    parms.lwd = o->lineWidth;
    parms.rbgcirc = 0.0;
    parms.spacing = o->markSpacing;
    mkgr_mark_grid_parms_round(&parms, pxpmm, 1.0);
    int32_t spacing_px = (int32_t)(parms.spacing*pxpmm + 0.5); /* Main mark spacing in pixels.*/
    fprintf(stderr, "pixels per cell = %d\n", spacing_px);
    
    mkgr_mark_grid_t *gr = mgs_create_mark_grid(HGX, HGY, &parms);
    
    /* Compute the image size: */
    int32_t NX = 2 * spacing_px * HGX;
    int32_t NY = 2 * spacing_px * HGY;
   
    /* Create the float image and fill it with black: */
    float_image_t *img = float_image_new(1,NX,NY);
    float_image_fill(img, 0.0f);
    
    r2_t ctr = (r2_t){{ 0.5*NX, 0.5*NY }}; /* Pixel coordinates of center of grid. */
    
    int32_t subsamp = 2; /* Antialiasing subsampling order. */
    mkgr_mark_grid_paint_image(img, 1, NULL, gr, pxpmm, &ctr, subsamp);
    mgs_write_png_image(o->outPrefix, img);
    mgs_write_mark_data(o->outPrefix, "png", gr);

    fprintf(stderr, "--- PNG image file written ---\n");

    float_image_free(img);
    mkgr_mark_grid_free(gr);
  }

mkgr_mark_grid_t *mgs_create_mark_grid
  ( int32_t HGX, 
    int32_t HGY,
    mkgr_mark_grid_parms_t *parms
  )
  { mkgr_mark_grid_t *gr = mkgr_mark_grid_new();
    
    auto mkgr_mark_t def_main_mark(int32_t ix, int32_t iy);
      /* Mark at cell corner point {(ix,iy)}: open cicle of radius {parms.rcircle}
        at corners of {5x5} cells, else upright cross of radius {parms.rcross}. */

    auto mkgr_mark_t def_redot_mark(int32_t ix, int32_t iy);
      /* Extra mark at cell corner point {(5*ix,5*iy)}: solid dot
        of radius {parms.rdot} inside the circle already there. */

    auto mkgr_mark_t def_sudot_mark(int32_t ix, int32_t iy);
      /* Extra mark at cell center {(ix+0.5,iy+0.5)},
        except near main circles: solid dot of radius {parms.rdot}. */
      
    frgb_t white = (frgb_t){{ 1.000f, 1.000f, 1.000f }};
    frgb_t black = (frgb_t){{ 0.000f, 0.000f, 0.000f }};
    
    double spacing = parms->spacing;
      
    if (parms->rbgcirc > 0)
      { r2_t ctr = (r2_t){{ 0.0, 0.0 }};
        mkgr_mark_t mk = mkgr_make_dot(ctr, parms->rbgcirc, black);
        mkgr_mark_grid_append_mark(gr, &mk);
      }

    /* Main marks - crosses every {spacing}, circles every {5*spacing}: */
    i2_t iMin1 = (i2_t){{ -HGX, -HGY }};
    i2_t iMax1 = (i2_t){{ +HGX, +HGY }};
    mkgr_mark_grid_append_marks(gr, iMin1, iMax1, def_main_mark);
      
    /* Secondary marks - extra dots inside main circles: */
    i2_t iMin2 = (i2_t){{ -HGX/5, -HGY/5 }};
    i2_t iMax2 = (i2_t){{ +HGX/5, +HGY/5 }};
    mkgr_mark_grid_append_marks(gr, iMin2, iMax2, def_redot_mark);

    /* Tertiary marks - dots inside each cell of size {spacing} */
    i2_t iMin3 = (i2_t){{ -HGX, -HGY }};
    i2_t iMax3 = (i2_t){{ +HGX-1, +HGY-1 }};
    mkgr_mark_grid_append_marks(gr, iMin3, iMax3, def_sudot_mark);

    return gr;
    
    /* Internal implementations: */
    
    mkgr_mark_t def_main_mark(int32_t ix, int32_t iy)
      { r2_t ctr = (r2_t){{ ix*spacing, iy*spacing }};
        bool_t cross = ((ix % 5) != 0) || ((iy % 5) != 0);
        if (cross)
          { return mkgr_make_cross(ctr, parms->rcross, 0.0, parms->lwd, white); }
        else
          { return mkgr_make_circle(ctr, parms->rcircle, parms->lwd, white); }
      }

    mkgr_mark_t def_redot_mark(int32_t ix, int32_t iy)
      { r2_t ctr = (r2_t){{ ix*5*spacing, iy*5*spacing }};
        return mkgr_make_dot(ctr, parms->rdot, white);
      }

    mkgr_mark_t def_sudot_mark(int32_t ix, int32_t iy)
      { r2_t ctr = (r2_t){{ (ix+0.5)*spacing, (iy+0.5)*spacing }}; /* ??? */
        double rad = parms->rdot; 
        int32_t rx = ((ix % 5) + 6) % 5; /* Get {ix+1} mod 5. */
        int32_t ry = ((iy % 5) + 6) % 5; /* Get {iy+1} mod 5. */
        bool_t omit = ((rx <= 1) && (ry <= 1));
        if (omit) { rad = 0.0; }
        return mkgr_make_dot(ctr, rad, white);
      }
  }

void mgs_write_mark_data(char *prefix, char *tag, mkgr_mark_grid_t *gr)
  { 
    char *fname = jsprintf("%s_%s.txt", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    mkgr_mark_grid_write(wr, gr);
    fclose(wr);
    free(fname);
  }

void mgs_get_grid_bounding_box(mkgr_mark_grid_t *gr, double mrg, r2_t *sMin, r2_t *sMax)
  {
    /* Include the backgroudn disk and add the margin: */
    for (uint32_t j = 0;  j < 2; j++)
      { sMin->c[j] = gr->pMin.c[j] - mrg;
        sMax->c[j] = gr->pMax.c[j] + mrg;
      }
  }

void mgs_write_png_image(char *prefix, float_image_t *img)
  { 
    char *fname = jsprintf("%s.png", prefix);
    double v0 = 0.0;
    double vM = 1.0;
    double gammaEnc = 1.0;
    double bias = 0.0;
    bool_t verbose = FALSE;
    float_image_write_gen_named(fname, img, image_file_format_PNG, v0, vM, gammaEnc, bias, verbose);
    free(fname);
  }

options_t *mgs_parse_options(int32_t argc, char **argv)
  {
    /* Initialize the argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    argparser_get_keyword(pp, "-gridSize");
    o->gridSize.c[0] = (int32_t)argparser_get_next_int(pp, 0, 2000);
    o->gridSize.c[1] = (int32_t)argparser_get_next_int(pp, 0, 2000);
    if (((o->gridSize.c[0] %2) != 0) || ((o->gridSize.c[1] %2) != 0))
      { argparser_error(pp, "grid width and height must be even."); }

    if (argparser_keyword_present(pp, "-markRadius"))
      { o->markRadius = epswr_round_dim_mm(argparser_get_next_double(pp, 0, 1000.0), o->printerDPI);  }
    else
      { o->markRadius = epswr_round_dim_mm(0.16, o->printerDPI); }

    if (argparser_keyword_present(pp, "-lineWidth"))
      { o->lineWidth = argparser_get_next_double(pp, 0, 1000.0);  }
    else
      { o->lineWidth = epswr_round_dim_mm(0.08, o->printerDPI); }

    if (argparser_keyword_present(pp, "-markSpacing"))
      { o->markSpacing = argparser_get_next_double(pp, 0, 1000.0); }
    else
      { o->markSpacing = epswr_round_dim_mm(1.0, o->printerDPI); }

    if (argparser_keyword_present(pp, "-printerDPI"))
      { o->printerDPI = argparser_get_next_double(pp, 36, 1000.0); }
    else
      { o->printerDPI = 0.00; }

    if (argparser_keyword_present(pp, "-imagePXPMM"))
      { o->imagePXPMM = argparser_get_next_double(pp, 1.0, 1000.0); }
    else
      { o->imagePXPMM = 0.00; }

    if (argparser_keyword_present(pp, "-outPrefix"))
      { o->outPrefix = argparser_get_next(pp); }
    else
      { o->outPrefix = NULL; }

    /* Get positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for leftover args: */
    argparser_finish(pp);

    return o;
  }

