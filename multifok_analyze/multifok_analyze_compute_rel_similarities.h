#ifndef multifok_analyze_compute_rel_similarities_H
#define multifok_analyze_compute_rel_similarities_H
/* Last edited on 2018-09-10 18:33:41 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <float_image.h>

void multifok_analyze_compute_rel_similarities
  ( int32_t NF, 
    int32_t ND, 
    double e2[], 
    double **s2P, 
    double **unP,
    bool_t debug
  );
  /* Computes relative similarity coefficients {s2[]} and fractional misfocusing
    indicators {un[]} from quadratic discrepancies {e2[]} (observed or fitted)
    for each of {ND} significant window positions and {NF} frames.
    
    Assumes that the vector {e2} has {NP=ND*NF} elements.  Allocates {s2} and {un} 
    with the same size and returns them in {*s2P} and {*unP}. Specifically,
    assumes that element with index {p = d*NF + f} of these vectors corresponds
    to significant window position {d} on the frame with index {f}.
    
    The procedure operates on each window position {d}.  Let {e2f[f]}, for {f} in {0..NF-1},
    be {e2[d*NF + f}; and similarly for {s2f[f],unf[f]}.
    
    The relative similarities {s2f[0..NF-1]} are the discrepancies {e2f[0..NF-1]},
    complemented and normalized so that the approximate maximum discrepancy 
    is mapped to 0 and the approximate minimum discrepancy is mapped to 
    1, with outliers clipped to those values.  
    
    The misfocus {unf[f]} is the difference
    between the frame index {f} of the observation and the (fractional) frame 
    index {fopt} such that {s2f[fops]} seems to be maximum among {s2f[0..NF-1]}.
    
    If {debug} is true, prints some debugging info. */


#endif
