/* See {multifok_analyze_extract_window.h}. */
/* Last edited on 2018-01-01 23:26:31 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <multifok_focus_op.h>

#include <multifok_analyze_extract_window.h>

double multifok_analyze_extract_window
  ( int32_t NW, 
    float_image_t *img, 
    int32_t dx, 
    int32_t dy,
    double noise,
    double fr[], 
    double w[]
  )
  {
    demand(img->sz[0] == 1, "image must be monochromatic");
    demand((NW % 2) == 1, "window size must be odd");

    /* Extract window samples: */
    int32_t NS = NW*NW;
    float va[NS];
    double rep = TRUE; /* Replicate border pixels (should not matter). */
    float_image_get_window_samples(img, 0, dx, dy, NW, NW, rep, va);
    
    /* Compute weighted sample average in window: */
    double sum_w = 1.0e-200;
    double sum_w_v = 0.0;
    for (int32_t s = 0; s < NS; s++)
      { double vs = (double)va[s];
        double ws = w[s];
        demand(ws > 0, "invalid weight");
        sum_w += ws;
        sum_w_v += ws*vs;
      }
    demand(! isnan(sum_w_v), "samples are {NAN}");
    double avg = sum_w_v/sum_w;
    
    /* Compute weighted sample deviation in window, accounting for noise: */
    double sum_w_d2 = 0.0;
    for (int32_t s = 0; s < NS; s++)
      { double ds = ((double)va[s]) - avg;
        double ws = w[s];
        sum_w_d2 += ws*(ds*ds + noise*noise);
      }
    double dev = sqrt(sum_w_d2/sum_w);
    assert(dev >= 0.99999*noise);
    
    /* Normalize samples: */
    for (int32_t s = 0; s < NS; s++) { fr[s] = ((va[s] - avg)/dev); }
    
    return dev;
  }
