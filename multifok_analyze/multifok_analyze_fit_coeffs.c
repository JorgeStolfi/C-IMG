/* See {multifok_analyze_fit_coeffs.h}. */
/* Last edited on 2024-12-21 13:59:23 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <lsq_array.h>

#include <multifok_analyze_fit_coeffs.h>

void multifok_analyze_fit_coeffs
  ( int32_t NP,
    int32_t NQ, 
    double terms[],
    double e2[],
    double **coeffP,
    double **m2P
  )
  {
    /* Compute the linear regression coefficients: */
    double *coeff = notnull(malloc(NQ*sizeof(double)), "no mem");
    bool_t verbose = FALSE;
    int32_t rank = lsq_array_fit(NP, NQ, 1, terms, e2, NULL, coeff, verbose);
    if (rank < NQ)
      { fprintf(stderr, "!! warning - least squares fit system has rank %d < %d\n", rank, NQ); }
    /* Allocates and computes the fitted values: */
    double *m2 = notnull(malloc(NP*sizeof(double)), "no mem");
    for (uint32_t p = 0;  p < NP; p++)
      { double sum = 0;
        for (uint32_t q = 0;  q < NQ; q++)
          { sum += coeff[q] * terms[p*NQ + q]; }
        m2[p] = sum;
      }
    
    (*coeffP) = coeff;
    (*m2P) = m2;
  }
