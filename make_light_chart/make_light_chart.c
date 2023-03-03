#define PROG_NAME "make_light_chart"
#define PROG_DESC "creates an EPS file with a simple reflectance chart"
#define PROG_VERS "1.0"

/* Last edited on 2023-03-02 22:42:46 by stolfi */

#define make_light_chart_C_COPYRIGHT \
  "Copyright © 2009  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    [ -rings {RINGS} ] \\\n" \
  "    [ -radius {RADIUS} ] \\\n" \
  "    [ -showSpotNumbers ] \\\n" \
  "    -outPrefix {OUT_PREFIX}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Writes and EPS file containing the a reflectance target, and a text description of the same.\n" \
  "\n" \
  "  The target consists of up to 25 circles with the given max radius" \
  " (in mm) against a black background.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -radius {RADIUS}\n" \
  "    Specifies the radius of each disk in millimeters. The disk itself\n" \
  " will be somewhat smaller.\n" \
  "\n" \
  "  -rings {RINGS}\n" \
  "    Specifies the number of concentric rings of spots (1 = only the central dot).\n" \
  "\n" \
  "  -showSpotNumbers\n" \
  "    If this optionis specified, the program will write the" \
  " spot index inside each spot of the chart.\n" \
  "\n" \
  "  -outPrefix {OUT_PREFIX}\n" \
  "    This mandatory parameter specifies the prefix for output filenames" \
  " minus the \".eps\" or \".txt\" extensions.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgmtopsx(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2009-07-03 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2021-01-02 J.SRolfi converted from {pswr.h} to {epswr.h}.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " make_light_chart_C_COPYRIGHT "\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <r2.h>
#include <jsfile.h>
#include <jsstring.h>
#include <epswr.h>
#include <epswr_iso.h>
#include <argparser.h>

typedef struct options_t
  { double radius;
    int32_t rings;
    char *outPrefix;
    bool_t showSpotNumbers;
  } options_t;

typedef struct chart_t 
  { double radius;         /* Radius of background-colored area. */
    double bg_lum;         /* Reflectivity of background area. */
    int32_t nspots;            /* Number of spots in chart. */
    r2_vec_t spot_ctr;     /* Spot centers. */
    double_vec_t spot_rad; /* Spot radii. */
    double_vec_t spot_lum; /* Spot reflectivities. */
  } chart_t;
  /* Geometry of the chart. */

/* PROTOTYPES */

int32_t main (int32_t argc, char **argv);
chart_t *compute_chart(int32_t rings, double d_spot[], double r_spot[], double r_targ, double rel_marg);
void draw_chart(epswr_figure_t *eps, chart_t *ch, bool_t showSpotNumbers);
void write_chart_description(FILE *wr, chart_t *ch);
options_t *parse_options(int32_t argc, char **argv);
  /* Obtains arguments from the command line. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { 
    /* Command-line parameters: */
    options_t *o = parse_options(argc, argv);
    
    int32_t rings = o->rings; /* Number of spot rings (including center spot). */
    
    /* Basic target geometry: */
    double r_shrink = 0.80;  /* Spot radius reduction from layer to layer. */
    double r_spot[rings];    /* Spot radius in each ring of spots (incl. margin). */
    double d_spot[rings];    /* Spot center distance from target center. */
    double c30 = cos(M_PI/6);
    double s30 = sin(M_PI/6);
    double c15 = cos(M_PI/12);
    double s15 = sin(M_PI/12);
    for (int32_t j = 0; j < rings; j++)
      { if (j == 0)
          { r_spot[j] = o->radius;
            d_spot[j] = 0.0;
          }
        else if (j == 1)
          { r_spot[j] = r_shrink*r_spot[j-1];
            d_spot[j] = d_spot[j-1] + r_spot[j-1] + r_spot[j];
          }
        else if (j == 2)
          { r_spot[j] = r_shrink*r_spot[j-1];
            double da = d_spot[j-2] + r_spot[j-2] + r_spot[j];
            double hb = d_spot[j-1]*s30;
            double cb = r_spot[j-1] + r_spot[j];
            double db = d_spot[j-1]*c30 + sqrt(fmax(0, cb*cb - hb*hb));
            d_spot[j] = fmax(da, db);
          }
        else if (j == 3)
          { r_spot[j] = r_shrink*r_spot[j-1];
            double ha = d_spot[j-2]*s15;
            double ca = r_spot[j-2] + r_spot[j];
            double da = d_spot[j-2]*c15 + sqrt(fmax(0, ca*ca - ha*ha));
            double hb = d_spot[j-1]*s15;
            double cb = r_spot[j-1] + r_spot[j];
            double db = d_spot[j-1]*c15 + sqrt(fmax(0, cb*cb - hb*hb));
            d_spot[j] = fmax(da, db);
          }
        fprintf(stderr, "ring %2d dist = %9.4f  radius = %9.4f\n", j, d_spot[j], r_spot[j]);
      }
    
    double r_targ = d_spot[rings-1] + r_spot[rings-1]; /* Rel. target radius (minus marg.) */
    double rel_marg = 1.10;  /* Relative margin to leave around each spot. */

    chart_t *ch = compute_chart(rings, d_spot, r_spot, r_targ, rel_marg);
    
    /* Open the figure stream/document: */
    double figRad_mm = rel_marg*r_targ;
    double figRad_pt = figRad_mm*epswr_pt_per_mm;
    double hSize = 2*figRad_pt;
    double vSize = 2*figRad_pt;
    double mrg = 4.0; /* Margin width (pt). */
    bool_t eps_verbose = FALSE;
    epswr_figure_t *eps = epswr_new_named_figure
      ( NULL, o->outPrefix, NULL, -1, NULL, 
        hSize, vSize, mrg, mrg, mrg, mrg,
        eps_verbose
      );
    double xMin = -figRad_mm, xMax = +figRad_mm;
    double yMin = -figRad_mm, yMax = +figRad_mm;
    epswr_set_client_window(eps, xMin, xMax, yMin, yMax);
    draw_chart(eps, ch, o->showSpotNumbers);
    epswr_end_figure(eps);

    /* Write chart geometry file: */
    char *geomFileName = txtcat(o->outPrefix, ".txt");
    FILE *geomFile = open_write(geomFileName, TRUE);
    write_chart_description(geomFile, ch);
    fclose(geomFile);
    free(geomFileName);
      
    return 0;
  }

chart_t *compute_chart(int32_t rings, double d_spot[], double r_spot[], double r_targ, double rel_marg)
  {
    chart_t *ch = (chart_t *)notnull(malloc(sizeof(chart_t)), "no mem");
    ch->bg_lum = 0.000;
    ch->radius = r_targ*rel_marg;
    ch->spot_ctr = r2_vec_new(25);
    ch->spot_rad = double_vec_new(25);
    ch->spot_lum = double_vec_new(25);
    int32_t max_index = -1;
    for (int32_t j = 0; j < rings; j++)
      { int32_t nsp;        /* Number of spots in ring. */
        int32_t k_base;     /* Lightness index of first spot. */
        int32_t k_step;     /* Increment in lightness index from spot to spot. */
        double t_base;  /* Angular position of first spot. */
        double t_step;  /* Angular increment between spots. */
        if (j == 0)
          { nsp = 1;  k_base = -1; k_step = 0; t_base = 0; t_step = 0; }
        else if (j == 1)
          { nsp = 6;  k_base = 0;  k_step = 4; t_base = 0; t_step = M_PI/3; }
        else if (j == 2)
          { nsp = 6;  k_base = 2;  k_step = 4; t_base = M_PI/6; t_step = M_PI/3; }
        else if (j == 3)
          { nsp = 12; k_base = 1;  k_step = 2; t_base = M_PI/12; t_step = M_PI/6; }
        else
          { assert(FALSE); }
        
        for (int32_t i = 0; i < nsp; i++)
          {
            double xc, yc;  /* Center of spot. */
            double rs;      /* Radius of spot. */
            double spY;     /* Lightness of spot. */
            /* Compute lightness index and lightness: */
            int32_t k = k_base + i*k_step;  /* Lightness index. */
            if (k == -1)
              { spY = 1 - ch->bg_lum; }
            else
              { spY = pow(0.5, (k+1)/9.0); }
            /* Compute actual position in circle: */
            int32_t it = (i < nsp/2 ? 2*i : 2*(i - nsp/2)+1);
            double t = t_base + it*t_step;
            double rc =  d_spot[j];
            xc = rc*cos(t);
            yc = rc*sin(t);
            rs = r_spot[j]/rel_marg;
            int32_t index = k+1;
            if (index > max_index) { max_index = index; }
            fprintf(stderr, "%2d %+9.4f %+9.4f %8.4f  %6.4f\n", index, xc, yc, rs, spY);
            r2_vec_expand(&(ch->spot_ctr), index);
            double_vec_expand(&(ch->spot_rad), index);
            double_vec_expand(&(ch->spot_lum), index);
	    ch->spot_ctr.e[index] = (r2_t){{ xc, yc }};
            ch->spot_rad.e[index] = rs;
            ch->spot_lum.e[index] = spY;
          }
      }
    ch->nspots = max_index + 1;
    r2_vec_trim(&(ch->spot_ctr), ch->nspots);
    double_vec_trim(&(ch->spot_rad), ch->nspots);
    double_vec_trim(&(ch->spot_lum), ch->nspots);
    return ch;
  }

void draw_chart(epswr_figure_t *eps, chart_t *ch, bool_t showSpotNumbers)
  { 
    epswr_set_pen(eps, 1.000, 0.000, 0.000,  0.40,  0.0, 0.0);
    epswr_set_label_font(eps, "Courier", 10.0);
    epswr_comment(eps, "background circle:");
    epswr_set_fill_color(eps, ch->bg_lum, ch->bg_lum, ch->bg_lum);
    epswr_circle(eps, 0.000, 0.000, ch->radius, TRUE, FALSE);
    for (int32_t i = 0; i < ch->nspots; i++)
      {
        char *cmt = NULL;
        asprintf(&cmt, "spot number %d:", i);
        epswr_comment(eps, cmt);
        free(cmt);
        r2_t ctr = ch->spot_ctr.e[i];
        double rad = ch->spot_rad.e[i];
        double spY = ch->spot_lum.e[i];

        epswr_set_fill_color(eps, spY, spY, spY);
        epswr_circle(eps, ctr.c[0], ctr.c[1], rad, TRUE, FALSE);
        if(showSpotNumbers)
          { double txY = (spY > 0.50 ? 0.00 : 1.00);
            epswr_set_fill_color(eps, txY, txY, txY);
            char* lab = NULL;
            asprintf(&lab, "%02d", i);
            epswr_label(eps, lab, "0", ctr.c[0],ctr.c[1],  0.0, TRUE, 0.5,0.5, TRUE, FALSE);
            free(lab);
          }
      }
  }

void write_chart_description(FILE *wr, chart_t *ch)
  {
    fprintf(wr,"nspots = %d\n",ch->nspots);
    fprintf(wr,"radius = %f\n",ch->radius);
    fprintf(wr,"bg_lum = %6.4f\n",ch->bg_lum);
    for (int32_t i = 0; i < ch->nspots; i++)
      { r2_t ctr = ch->spot_ctr.e[i];
        double rad = ch->spot_rad.e[i];
        double spY = ch->spot_lum.e[i];
        fprintf(wr, "%2d %+9.4f %+9.4f %8.4f  %6.4f\n", i, ctr.c[0], ctr.c[1], rad, spY);
      }
    fflush(wr);
  }

options_t *parse_options(int32_t argc, char **argv)
  {
    /* Initialize the argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    if (argparser_keyword_present(pp, "-radius"))
      { o->radius = argparser_get_next_double(pp, 0, 1000.0);  }
    else
      { o->radius = 3.0; }

    if (argparser_keyword_present(pp, "-rings"))
      { o->rings = (int32_t)argparser_get_next_int(pp, 1, 4);  }
    else
      { o->rings = 4; }

    o->showSpotNumbers = argparser_keyword_present(pp, "-showSpotNumbers");

    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);

    /* Get positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for leftover args: */
    argparser_finish(pp);

    return o;
  }
