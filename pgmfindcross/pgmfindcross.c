#define PROG_NAME "pgmfindcross"
#define PROG_DESC "sub-pixel location of crosshair in PGM image"
#define PROG_VERS "1.0"

/* Last edited on 2023-02-25 16:05:02 by stolfi */
/* Copyright © 2003 by the State University of Campinas (UNICAMP).*/
/* See the copyright, authorship, and warranty notice at end of file. */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [ -linewd LW ] [ -radius R ] \\\n" \
  "  [ -maxdisp MD ] [ -maxrot MW ] \\\n" \
  "  PGMFILE"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "???\n"\
  "OPTIONS\n" \
  "???\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "???\n" \
  "\n" \
  "AUTHOR\n" \
  "???"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <values.h>

#include <affirm.h>
#include <argparser.h>
#include <jspnm.h>
#include <jsmath.h>
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>
#include <rn.h>

/* Number of {h,v} steps per pixel: */
#define NSTEPS (10)

typedef struct options_t
  { char* infile;     /* Image file name. */
    int radius;       /* Window radius. */
    double linewd;    /* Nominal line width, in pixels. */
    double maxdisp;   /* Maximum coordinate adjustment (pixels). */
    double maxrot;    /* Maximum direction adjustment (degrees). */
    bool_t darklines; /* TRUE if lines are darker than background. */
  } options_t;
  
typedef struct Tables
  { int prad;     /* The profile {Pr(dp/NSTEPS)} is zero for {dp > prad}. */
    double *pr;   /* Line profile: {pr[abs(dp)] = Pr(dp/NSTEPS)}. */
    int wrad;     /* The weight {Wt(hd/NSTEPS)} is zero for {hd > wrad}. */
    double *wt;   /* Projected weight: {wt[abs(wt)] = Wt(hd/NSTEPS)}. */
  } Tables;

extern double strtod(const char *, char **);
int main(int argc, char **argv);
options_t *parse_options(int argc, char **argv);
 
void find_crossing
  ( uint16_image_t *img,
    double *hp, double *vp,  /* (IN/OUT) Coordinates of crossing (pixels). */
    double *at, double *bt,  /* (IN/OUT) Directions of lines (degrees). */
    options_t *o,
    Tables *tb
  );
  /* Given the approximate coordinates {hp,vp} of a line-crossing,
    and approximate directions {at,bt} of the lines, finds
    more precise values {*hp,*vp,*at,*bt} for those
    parameters. */
    
void find_line
  ( uint16_image_t *img,
    int *hp, int *vp,  /* (IN/OUT) Coordinates of point on line, times {NSTEPS}. */
    double *t,         /* (IN/OUT) Direction of line (degrees). */
    int dp_max,        /* Number of alt. positions to try on either side. */
    int p_step,        /* Position step, in {1/NSTEP} pixels. */
    int dt_max,        /* Number of alt.directions to try on either side. */
    double t_step,     /* Direction step, in degrees. */
    options_t *o,
    Tables *tb
  );
  /* Given the approximate coordinates of a point {*hp,*vp} on a
    line, and its approximate direction {*t} refines those
    parameters by tilting the line and displacing the point in the
    transversal direction. The angle perturbation ranges from {-dt_max}
    to {+dt_max} times {t_step} degrees; the position perturbation
    ranges from {-dp_max} to {+dp_max} times {p_step/NSTEPS} pixels. */

double match_line
  ( uint16_image_t *img,
    int hp, int vp,       /* Coordinates of window center, times {NSTEPS}. */
    double hn, double vn, /* Unit vector normal to line. */
    options_t *o,
    Tables *tb
  );
  /* Computes the dot product of the image and the ideal line that
    goes through the point {(hp/NSTEPS,vp/NSTEPS)} and is orthogonal
    to the vector {(hn,vn)}. */
 
double ideal_pixel(int hd, int vd, double hn, double vn, options_t *o, Tables *tb);
  /* Computes the ideal value of a pixel whose center lies
   {(hd/NSTEPS,vd/NSTEPS)} away from the axis of a line 
   orthogonal to the vector {(hn,vn)}, whose
   transversal profile is given by {tb->pr}. */
   
double window_weight(int hd, int vd, options_t *o, Tables *tb);
  /* Computes the window weight function for a pixel whose 
   center lies {(hd/NSTEPS,vd/NSTEPS)} away from the window center. */
 
Tables *make_tables(options_t *o, int maxval);
  /* Computes auxiliary tables (window weights, line profile, etc.) */

double *make_half_gaussian
  ( double sigma, 
    double eps, 
    int *radP
  );
  /* Tabulates a Gaussian-like weight function with deviation {sigma}.
    The weights are modified so that they become zero at when the true
    Gaussian would be less than {eps}. The profile is returned as an
    array {g[0..rad]}, sampled every {1/NSTEPS} pixels; where {rad =
    *radP} is a large enough integer, computed by the procedure. */

int main(int argc, char **argv)
  {
    options_t *o = parse_options(argc, argv);       /* Command-line options. */
    uint16_image_t *img = uint16_image_read_pnm_named(o->infile, FALSE);
    Tables *tb = make_tables(o, img->maxval);
    double hp, vp, at, bt;
    int nlin = 0, res;
    
    while ((res = scanf("%lg %lg %lg %lg", &hp, &vp, &at, &bt)) != EOF)
      { nlin++;
        if (res != 4) { fprintf(stderr, "stdin:%d: format error\n", nlin); exit(1); }
        find_crossing(img, &hp,&vp, &at,&bt, o, tb);
        printf("%6.1f %6.1f  %6.1f %6.1f\n", hp, vp, at, bt);
        fflush(stdout);
      }
    fclose(stdout);
    return 0;
  }

static int debugit = 0;

/* Minimum step in line direction (degrees): */
#define MINTSTEP (1.0)

/* Minimum step in crossing position ({1/NSTEPS} of pixel): */
#define MINPSTEP (1)

void find_crossing
  ( uint16_image_t *img,
    double *hp, double *vp,  /* (IN/OUT) Coordinates of crossing (pixels). */
    double *at, double *bt,  /* (IN/OUT) Directions of lines (degrees). */
    options_t *o,
    Tables *tb
  )
  { /* Convert coordinates into multiples of {1/NSTEPS} pixels: */
    int hpi = (int)((*hp)*NSTEPS + 0.5);
    int vpi = (int)((*vp)*NSTEPS + 0.5);
    /* Initial steps in position and direction: */
    int p_step = (o->linewd < 1.0 ? NSTEPS : (int)(ceil(o->linewd*NSTEPS)));
    double t_step = 4.0;
    /* Initial counts of position and direction trials: */
    int dp_max = (int)(floor((o->maxdisp*NSTEPS)/p_step));
    int dt_max = (int)(floor(o->maxrot/t_step));
    fprintf(stderr, 
      "--- begin find_cross (%6.1f, %6.1f,  %6.1f, %6.1f) ---\n",
      *hp, *vp, *at, *bt
    );
    /* Multiscale search (should vary sampling resolution...) */
    do
      { /* Current best positions of lines: */
        int hps, vps;
        if ((o->maxdisp*NSTEPS < MINPSTEP) || (p_step < MINPSTEP)) { dp_max = 0; }
        if ((o->maxrot < MINTSTEP) || (t_step < MINTSTEP)) { dt_max = 0; }
        fprintf(stderr, 
          "  current pos = (%4d/%d, %4d/%d) dir = (%6.1f, %6.1f)\n",
          hpi,NSTEPS, vpi,NSTEPS, *at, *bt
        );
        fprintf(stderr, 
          "  pos steps = (±%2d * %4d/%d) dir steps = (±%2d * %7.2f)\n",
          dp_max, p_step,NSTEPS, dt_max, t_step
        );
    
        /* Search positions of lines, independently: */
        hps = hpi; vps = vpi;
        fprintf(stderr, "  searching for first line...\n");
        find_line(img, &hps, &vps, at, dp_max, p_step, dt_max, t_step, o, tb);
        fprintf(stderr, "  searching for second line...\n");
        find_line(img, &hps, &vps, bt, dp_max, p_step, dt_max, t_step, o, tb);
        hpi = hps; vpi = vps;
        /* Next scale of trials: */
        p_step /= 2;
        t_step /= 2.0;
        dp_max = 2;
        dt_max = 2;
      }
    while ((p_step >= MINPSTEP) || (t_step >= MINTSTEP));
    (*hp) = ((double)hpi)/NSTEPS;
    (*vp) = ((double)vpi)/NSTEPS;
    fprintf(stderr, 
      "--- end find_cross (%6.1f, %6.1f,  %6.1f, %6.1f) ---\n", 
      *hp, *vp, *at, *bt
    );
  } 

void find_line
  ( uint16_image_t *img,
    int *hp, int *vp,  /* (IN/OUT) Coordinates of point on line, times {NSTEPS}. */
    double *t,         /* (IN/OUT) Direction of line (degrees). */
    int dp_max,        /* Number of alt. positions to try on either side. */
    int p_step,        /* Position step, in {1/NSTEP} pixels. */
    int dt_max,        /* Number of alt.directions to try on either side. */
    double t_step,     /* Direction step, in degrees. */
    options_t *o,
    Tables *tb
  ) 
  { double hp_in = (double)(*hp);
    double vp_in = (double)(*vp);
    double t_in = (*t);
    double best_score = -INF;
    int kp, kt, sp, st;
    st = 1;
    for (kt = 0; kt <= 2*dt_max; kt++) 
      { /* Signed direction perturbation: */
        int dt = st*((kt+1)/2); 
        /* Candidate inclination of line: */
        double t_c = t_in + t_step * (double)dt;
        /* Unit vector perpendicular to line: */
        double hn_c = sin((t_c)*M_PI/180.0);
        double vn_c = cos((t_c)*M_PI/180.0);
        sp = 1;
        for (kp = 0; kp <= 2*dp_max; kp++)
          { /* Signed position perturbation: */
            int dp = sp*((kp+1)/2); 
            /* Candidate line displacement in pixels: */
            double disp_c = p_step * (double)dp;
            /* Candidate center of line segment: */
            int hp_c = (int)(hp_in + disp_c * hn_c + 0.5);
            int vp_c = (int)(vp_in + disp_c * vn_c + 0.5);

            /* int trash = (debugit = ((dt == 0) & (dp == 0))); */
            double score = match_line(img, hp_c, vp_c, hn_c, vn_c, o, tb);
            fprintf(stderr, "  p = (%4d/%d %4d/%d ) t = %6.1f  score = %8.3f",
              hp_c,NSTEPS, vp_c,NSTEPS, t_c, score
            );
            if (score > best_score) 
              { best_score = score; 
                (*hp) = hp_c; (*vp) = vp_c; (*t) = t_c;
                fprintf(stderr, " (improved)");
              }
            fprintf(stderr, "\n");
            sp = -sp;
          }
        st = -st;
      }
  }

double match_line
  ( uint16_image_t *img,
    int hp, int vp,       /* Coordinates of window center, times {NSTEPS}. */
    double hn, double vn, /* Unit vector normal to line. */
    options_t *o,
    Tables *tb
  )
  { int rows = img->rows, cols = img->cols;
    assert(img->chns == 1); /* For now... */
    int HN = NSTEPS/2;
    
    /* Indices of image pixels covered by window: */
    int hi_min = (hp - tb->wrad + HN + 1)/NSTEPS;
    int hi_max = (hp + tb->wrad + HN)/NSTEPS;
    int vi_min = (vp - tb->wrad + HN + 1)/NSTEPS;
    int vi_max = (vp + tb->wrad + HN)/NSTEPS;
    
    int hi, vi;
    double mv = (double)(img->maxval);
    double swyp = 0.0, swpp = 0.0, swy = 0.0, swp = 0.0, sw = 0.0;
    double dot;

    assert(NSTEPS == 2*HN);
    if (debugit) 
      {  fprintf(stderr, 
           "    --- begin match_mask(%4d/%d, %4d/%d,  %7.4f, %7.4f) ---\n", 
           hp,NSTEPS, vp,NSTEPS, hn, vn
         );
      }
      
    if (hi_min < 0) { hi_min = 0; }
    if (hi_max >= cols) { hi_max = cols-1; }
    if (vi_min < 0) { vi_min = 0; }
    if (vi_max >= rows) { vi_max = rows-1; }
    
    for (vi = vi_min; vi <= vi_max; vi++)
      { uint16_t *pixv = img->smp[vi];
        /* V displacement from window center {hp,vp} to pixel center: */
        int vd = vi*NSTEPS + HN - vp;
        if (debugit && (abs(vd) < 5*NSTEPS)) 
          { fprintf (stderr, "     "); }
        for (hi = hi_min; hi <= hi_max; hi++)
          { int hd = hi*NSTEPS + HN - hp;
            double w = window_weight(hd, vd, o, tb);
            double y;
            if (w > 0.0)
              { double p = ideal_pixel(hd, vd, hn, vn, o, tb);
                y = ((double)(pixv[hi]))/mv;
                swyp += w*y*p;
                swpp += w*p*p;
                swy += w*y;
                swp += w*p;
                sw += w;
              }
            else
              { y = 0.0; }
            if (debugit && (abs(vd) < 5*NSTEPS) && (abs(hd) < 5*NSTEPS))
              { if (w > 0) 
                  { fprintf (stderr, " %5.2f", y); }
                else
                  { fprintf (stderr, " %5s", ""); }
              }
          }
        if (debugit && (abs(vd) < 5*NSTEPS)) { fprintf (stderr, "\n"); }
      }
    /* Normalize ideal pixels {p} to zero mean: */
    if (sw != 0.0) 
      { double mp = swp/sw;
        swp = 0.0;
        swyp = swyp - mp*swy;
        swpp = swpp - 2.0*mp*swp + mp*mp*sw;
        if (swpp < 0.0) { swpp = 0.0; swp = 0.0; swyp = 0.0; }
      }
    /* Normalize ideal pixels to unit norm: */
    if (swpp != 0.0)
      { double zp = sqrt(swpp/sw);
        swp /= zp;
        swyp /= zp;
        swpp = sw;
      }
    dot = swyp/sw;
    /*
    fprintf (stderr, "    sw = %5.3f  avg: yp = %5.3f", sw, swyp/sw);
    fprintf (stderr, " p = %5.3f y = %5.3f\n", swp/sw, swy/sw);
    */
    if (o->darklines) { dot = -dot; }
    if (debugit) 
      { fprintf(stderr, 
          "    --- end match_mask(%4d/%d, %4d/%d,  %7.4f, %7.4f) dot = %6.4f ---\n",
          hp,NSTEPS, vp,NSTEPS, hn, vn, dot
        );
      }
    debugit = 0;
    return dot;
  }

Tables *make_tables(options_t *o, int maxval)
  { Tables *tb = (Tables *)pnm_malloc(sizeof(Tables));
    tb->pr = make_half_gaussian(o->linewd/2.0, 0.005, &(tb->prad));
    fprintf(stderr, "profile radius = %d/%d\n", tb->prad, NSTEPS);
    tb->wt = make_half_gaussian(o->radius, 0.02, &(tb->wrad));
    fprintf(stderr, "window radius = %d/%d\n", tb->wrad, NSTEPS);
    return tb;
  }

double window_weight(int hd, int vd, options_t *o, Tables *tb)
  { if (hd < 0) { hd = -hd; }
    if (vd < 0) { vd = -vd; }
    if (hd*hd + vd*vd > tb->wrad*tb->wrad)
      { return 0.0; }
    else
      { return tb->wt[hd]*tb->wt[vd]; }
  }

double ideal_pixel(int hd, int vd, double hn, double vn, options_t *o, Tables *tb)
  { double t = fabs(((double)hd)*hn + ((double)vd)*vn);
    int tk = (int)(floor(t));
    if (tk > tb->prad)
      { return 0.0; }
    else 
      { double p0 = tb->pr[tk];
        double tr = t - (double)tk;
        if (tr == 0.0)
          { return p0; }
        else 
          { double p1 = (tk >= tb->prad ? 0.0 : tb->pr[tk+1]);
            assert((tr > 0.0) && (tr < 1.0));
            return (1.0-tr)*p0 + tr*p1;
          }
      }
  }
  
double *make_half_gaussian
  ( double sigma, 
    double eps, 
    int *radP
  )
  { /* Distance at which the profile is equal to {eps}: */
    double maxdist = sqrt(-2.0*log(eps))*sigma;
    /* Number of profile samples needed to cover in {maxdist}: */
    int rad = (int)(ceil(maxdist*NSTEPS));
    /* The profile vector: */
    double *g = rn_alloc(rad + 1);
    int i;
    assert(g != NULL);
    /* Compute profile (a gaussian with deviation  {linewd/2.0}): */
    for (i = 0; i <= rad; i++)
      { double d = ((double)i)/NSTEPS;
        double t = d/sigma;
        double gauss = exp(-t*t/2.0);
        double h = (d/maxdist);
        double p = 1.0 - h*h*h*h;
        g[i] = p*p*gauss;
      }
    (*radP) = rad;
    return g;
  }

options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
            
    o->darklines = argparser_keyword_present(pp, "-darklines");

    if (argparser_keyword_present(pp, "-linewd"))
      { o->linewd = argparser_get_next_double(pp, 1.0e-100, 100.0); }
    else
      { o->linewd = 1.0; }
    
    if (argparser_keyword_present(pp, "-radius"))
      { o->radius = (int)ceil(argparser_get_next_double(pp, 0.0, 1000.0)); }
    else
      { o->radius = (int)(ceil(3.0*o->linewd) + 2); }
    
    if (argparser_keyword_present(pp, "-maxdisp"))
      { o->maxdisp = argparser_get_next_double(pp, 1.0e-6, 100.0); }
    else
      { o->maxdisp = 10.0; }
    
    if (argparser_keyword_present(pp, "-maxrot"))
      { o->maxrot = argparser_get_next_double(pp, 0.0, 90.0); }
    else
      { o->maxrot = 10.0; }

    argparser_skip_parsed(pp);
    
    if (argparser_next(pp) != NULL)
      { o->infile = argparser_get_next(pp); }
    else
      { o->infile = "-"; }

    argparser_finish(pp);
        
    return o;
  }

/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 2003 by the State University of Campinas (UNICAMP).
**
** Created feb/2003 by Amandio Sena Jr., IC-UNICAMP.
** Modified 2003-2004 by Jorge Stolfi, IC-UNICAMP.       
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
