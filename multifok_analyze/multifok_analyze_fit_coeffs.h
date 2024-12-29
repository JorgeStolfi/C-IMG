#ifndef multifok_analyze_fit_coeffs_H
#define multifok_analyze_fit_coeffs_H
/* Last edited on 2024-12-21 13:59:19 by stolfi */

#include <stdint.h>

void multifok_analyze_fit_coeffs
  ( int32_t NP,
    int32_t NQ, 
    double terms[],
    double e2[],
    double **coeffP,
    double **m2P
  );
  /* Assumes that {terms} is a matrix of quadratic terms, with {NP}
    rows and {NQ} columns, linearized by rows, and {e2} is a column
    vector of quadratic discrepancies, with {NP} elements.
    
    Each row {p} of {terms} is assumed to be observed values of {NQ}
    independent variables {Q[p][0..NQ-1]}, and the
    element {e2[p]} is assumed to the the corresponding 
    observed value of the dependendt variable to be fitted.

    Allocates and computes a vector {coeff[0..NQ-1]} of coefficients of the linear combination
    of columns of {terms} that best explains the quadratic discrepancies
    {e2[0..NP-1]}.  Returns {coeff} in {*coeffP}.
    
    Also allocates and computes a vector of fitted values {m2[0..NP-1]}, one for
    each observation.  Namely. {m2[p] = SUM{coeff[q]*Q[p][q]}}. Returns
    the vector {m2} in {*m2P}. */ 

#endif
