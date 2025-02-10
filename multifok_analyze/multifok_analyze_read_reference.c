/* See {multifok_analyze_read_reference.h}. */
/* Last edited on 2025-01-30 07:39:20 by stolfi */

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

#include <multifok_analyze_read_reference.h>

float_image_t *multifok_analyze_read_reference
  ( char *fname, 
    image_file_format_t ffmt, 
    bool_t verbose
  )
  {
    bool_t yUp = TRUE;
    float v0 = 0.0;
    float vM = 1.0;
    double gammaDec, bias; /* Enconding gammaDec and bias specified or implied by file. */
    float_image_t *fimg = float_image_read_gen_named(fname, ffmt, yUp, v0, vM, NULL, &gammaDec, &bias, verbose);
    demand(fimg->sz[0] == 1, "reference image should be monochromatic");
    return fimg;

  }
