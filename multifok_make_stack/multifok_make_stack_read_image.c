/* See {multifok_make_stack_read_image.h}. */
/* Last edited on 2025-01-30 05:04:23 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <vec.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_read_gen.h>
#include <image_file_format.h>

#include <multifok_make_stack_read_image.h>

float_image_t *multifok_make_stack_read_image
  ( char *fname, 
    image_file_format_t ffmt, 
    bool_t yUp, 
    bool_t verbose
  )
  {
    float v0 = 0.0;
    float vM = 1.0;
    double gammaDec, bias; /* Enconding gammaDec and bias specified or implied by file. */
    float_image_t *fimg = float_image_read_gen_named(fname, ffmt, yUp, v0, vM, NULL, &gammaDec, &bias, verbose);
    demand(fimg->sz[0] == 1, "reference image should be monochromatic");
    return fimg;
  }
