#ifndef multifok_analyze_extract_window_H
#define multifok_analyze_extract_window_H
/* Last edited on 2017-12-27 16:24:31 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <float_image.h>

double multifok_analyze_extract_window
  ( int32_t NW, 
    float_image_t *img, 
    int32_t dx, 
    int32_t dy,
    double noise,
    double fr[],
    double w[]
  );
  /* Assumes that {img} is a monochromatic image 
  
    Extracts the {NS=NW*NW} samples of the {NW} by {NW} window centered at the
    pixel in column {dx} and row {dy} of {img}.
    
    Then the procedure normalizes those {NS} values to zero mean and
    unit deviation, with pixel weights {w[0..NS-1]}. The deviation is
    biased by {noise}. In particular, if the image is constant in the
    window, the deviation is assumed to be {noise}, and the
    normalization maps all samples to zero.
    
    The result is returned in {fr[0..NS-1]}.
    
    The return value is the weighted standard deviation
    of the samples, before normalization and remapping,
    biased by {noise}. */
  
#endif

