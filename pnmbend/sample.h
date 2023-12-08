/* Last edited on 2017-06-22 17:59:58 by stolfilocal */

/* Subpixel sampling of a Gaussian-filtered function. */

#ifndef sample_H
#define sample_H

#define FILTER_RADIUS (4)
#define FILTER_SIGMA (2.0)

void sample_pixel(float_image_buffer_t *im, double xp, double yp, int chan, double *pp);
  /* Stores in `*pp' the value of channel `chan' of image `im' computed
    at the point with row `yp' and column `xp´ (which may be fractional).
    Assumes the band of the image contained in buffer `im' extends
    at least `FILTER_RADIUS' rows above and below of that point (if those rows
    exist).
  */

#endif
