/* See {multifok_analyze_compute_rel_similarities.h}. */
/* Last edited on 2018-09-10 18:38:00 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <wt_table.h>

#include <multifok_analyze_compute_rel_similarities.h>

void multifok_analyze_compute_rel_simil_aux
  ( int32_t NF, 
    double e2f[], 
    double s2f[], 
    double unf[],
    bool_t debug
  );
  /* Computes relative similarity coefficients {s2f[0..NF-1]} and fractional misfocusing
    indicators {un[0..NF-1]} for one significant window position,
    from the given quadratic discrepancies {e2f[0..NF-1]}.
    
    The relative similarities {s2f[0..NF-1]} are the discrepancies {e2f[0..NF-1]},
    complemented and normalized so that the approximate maximum discrepancy 
    is mapped to 0 and the approximate minimum discrepancy is mapped to 
    1, with outliers clipped to those values.  
    
    The misfocus {unf[f]} is the difference
    between the frame index {f} of the observation and the (fractional) frame 
    index {fopt} such that {s2f[fops]} seems to be maximum among {s2f[0..NF-1]}.
    
    If {debug} is true, prints some debugging info. */

void multifok_analyze_compute_rel_similarities
  ( int32_t NF, 
    int32_t ND, 
    double e2[], 
    double **s2P, 
    double **unP,
    bool_t debug
  )
  {
    int32_t NP = ND*NF; /* Number of observations. */
    
    double *s2 = notnull(malloc(NP*sizeof(double)), "no mem");    /* The relative similarities, with {NP} elements. */
    double *un = notnull(malloc(NP*sizeof(double)), "no mem");   /* The estimated misfocus, with {NP} elements. */
    for (int32_t d = 0; d < ND; d++)
      { double *e2f = &(e2[d*NF]); 
        double *s2f = &(s2[d*NF]); 
        double *unf = &(un[d*NF]); 
        bool_t debug_rel_sim = debug && ((d % 100) == 17);
        multifok_analyze_compute_rel_simil_aux(NF, e2f, s2f, unf, debug_rel_sim);
      }
    (*s2P) = s2;
    (*unP) = un;
  }


void multifok_analyze_compute_rel_simil_aux
  ( int32_t NF, 
    double e2f[], 
    double s2f[], 
    double unf[],
    bool_t debug
  )
  { /* Find maximum and minimum of {e2f[0..NF-1]}: */
    double e2max = -INF;
    double e2min = +INF;
    for (int32_t f = 0; f < NF; f++)
      { if (e2f[f] > e2max) { e2max = e2f[f]; }
        if (e2f[f] < e2min) { e2min = e2f[f]; }
      }
    if (debug) { fprintf(stderr, "    discrepancy  min %18.9f  max %18.9f\n", e2min, e2max); }
    
    double e2inf = e2min; /* Bottom value to complement to. */
    double e2sup = e2max; /* Top value to complement to. */

    /* Main iteration: */
    int32_t niter = 3;
    double favg = NAN;
    for (int32_t iter = 0; iter < niter; iter++)
      { 
        if (debug) { fprintf(stderr, "    --- iteration %d ---\n", iter); }
      
        /* Complement the discrepancies: */
        if (debug) { fprintf(stderr, "    complementing relative to [ %18.9f _ %18.9f ]\n", e2inf, e2sup); }
        for (int32_t f = 0; f < NF; f++) 
          { s2f[f] = (e2sup - e2f[f])/(e2sup - e2inf);
            if (s2f[f] < 0.0) { s2f[f] = 0.0; }
            if (s2f[f] > 1.0) { s2f[f] = 1.0; }
          }

        /* Find average index and its standard deviation: */
        favg = wt_table_avg(NF, s2f);
        assert(! isnan(favg));
        double fdev = sqrt(wt_table_var(NF, s2f, favg));
        assert(! isnan(fdev));
        if (debug) { fprintf(stderr, "    peak at %8.3f radius %8.3f\n", favg, fdev); }
        
        /* Estimate the width {frad} of the dip: */
        double frad = 3*fdev;
        if (frad > 10.0) { frad = 10.0; }
        if (frad < 2.0) { frad = 2.0; }
        if (debug) { fprintf(stderr, "    assumed radius %8.3f\n", frad); }

        /* Compute the average {e2avg} of {e2f[f]} ignoring {favg ± frad}: */
        double sum_we2 = 0.0;
        double sum_w = 1.0e-200;
        for (int32_t f = 0; f < NF; f++) 
          { if (fabs(((double)f) - favg) > frad) 
              { sum_we2 += e2f[f];
                sum_w += 1.0;
              }
          }
        double e2avg = sum_we2/sum_w;

        /* Compute the deviation {e2dev} about {e2avg} ignoring {favg ± frad}: */
        double sum_wd2 = 0.0;
        for (int32_t f = 0; f < NF; f++) 
          { if (fabs(((double)f) - favg) > frad) 
              { double d = e2f[f] - e2avg;
                sum_wd2 += d*d;
              }
          }
        double e2dev = sqrt(sum_wd2/sum_w);
        if (debug) { fprintf(stderr, "    off-peak sqr discrepancy  avg %12.6f  dev %12.6f\n", e2avg, e2dev); }

        /* Redefine the top range {e2sup} as {e2avg} minus multiple of {e2dev}: */
        e2sup = e2avg - 2 * e2dev;
        if (debug) { fprintf(stderr, "    raw reference level %18.9f\n", e2sup); }
        /* But don't let it get close or below {e2inf}: */
        double e2alt = e2inf + 0.75*(e2max - e2inf);
        if (e2sup < e2alt) { e2sup = e2alt; }
        if (debug) { fprintf(stderr, "    adj reference level %18.9f\n", e2sup); }
      }
    if (debug) { fprintf(stderr, "\n"); }

    /* Estimate the misfocusings, assuming frame {favg} is the best focused one: */
    for (int32_t f = 0; f < NF; f++) { unf[f] = ((double)f) - favg; }
  }
    
