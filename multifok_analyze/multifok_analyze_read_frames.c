/* See {multifok_analyze_read_frames.h}. */
/* Last edited on 2018-01-03 14:21:27 by stolfilocal */

#define _GNU_SOURCE
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

#include <multifok_analyze_read_frames.h>

float_image_t **multifok_analyze_read_frames
  ( char *framePattern, 
    int32_t NF, 
    int32_t frameID[], 
    image_file_format_t ffmt, 
    int32_t NC,
    int32_t NX, 
    int32_t NY,
    bool_t verbose
  )
  {
    float_image_t **fimg = notnull(malloc(NF*sizeof(float_image_t*)), "no mem");
    
    for (int32_t f = 0; f < NF; f++)
      { char *fname = NULL;
        asprintf(&fname, framePattern, frameID[f]);
        if (verbose) { fprintf(stderr, "reading frame[%3d]  from file \"%s\"\n", f, fname); }
        float v0 = 0.0;
        float vM = 1.0;
        double gammaDec, bias; /* Enconding gammaDec and bias specified or implied by file. */
        float_image_t *fimgf = float_image_read_gen_named(fname, ffmt, v0, vM, NULL, &gammaDec, &bias, verbose);
        demand (fimgf->sz[0] == NC, "inconsistent frame channel count");
        demand (fimgf->sz[1] == NX, "inconsistent frame width");
        demand (fimgf->sz[2] == NY, "inconsistent frame height");
        if (verbose) { fprintf(stderr, "read %s  NC = %d NX = %d NY = %d ...\n", fname, NC, NX, NY); }
        demand(NC == 1, "frames should be monochromatic");
        fimg[f] = fimgf;
        free(fname);
      }
    return fimg;
  }

    
