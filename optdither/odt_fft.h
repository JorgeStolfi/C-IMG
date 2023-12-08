/* Fourier tools for optimum dither seeking. */
/* Last edited on 2023-03-18 09:42:43 by stolfi */

#ifndef odt_fft_H
#define odt_fft_H

#define _GNU_SOURCE
#include <stdint.h>

#include <float_image.h>

void odt_fft_filter(float_image_t *xf, double wtx[], double wty[]);
  /* Assumes {xf} is the Hartley transform of some image. Reduces
    the low-frequency components by a high-pass weight with
    complemented axis transfer functions {wtx} and {wty}. */

void odt_fft_permutize(float_image_t *xf);
  /* Converts the grayscale image {xf} to a (float) dither matrix
    by replacing every sample {s} with {(r(s) + 0.5)/ns}
    where {r(s)} is the rank of {s} in increasing order of value, 
    and {ns} is the number of pixels in {xf}. */

#endif
