/* See fvw_texture.h */
/* Last edited on 2025-01-21 10:37:20 by stolfi */

#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <GL/glu.h>

#include <float_image.h>
#include <frgb.h>
#include <frgb_path.h>
#include <affirm.h>

#include <fvw_texture.h>

frgb_t fvw_texture_color_from_value(float v, float cm_vmin, float cm_vmax)
  {
    /* Determine the ranges for the {v} to {z} mapping: */
    demand(cm_vmin < cm_vmax, "invalid value range");
    double z_min, z_max;
    if ((cm_vmin < 0) && (cm_vmax > 0))
      { float cm_vex = fmaxf(fabsf(cm_vmin), fabsf(cm_vmax));
        cm_vmin = -cm_vex;
        cm_vmax = +cm_vex;
        z_min = -1; z_max = +1;
      }
    else if (cm_vmin < 0)
      { z_min = -1; z_max = 0; }
    else if (cm_vmax > 0)
      { z_min = 0; z_max = +1; }
    
    /* Compute color index {z} in {[z_min _ z_max]}: */
    double z;
    if (v < cm_vmin)
      { z = z_min; }
    else if (v > cm_vmax)
      { z = z_max; }
    else if (cm_vmin == cm_vmax)
      { z = (z_min + z_max)/2; }
    else 
      { z = z_min + (z_max - z_min)*(v - cm_vmin)/(cm_vmax - cm_vmin); }
      
    /* Now convert {z} to an RGB color: */
    frgb_t clr;
    if (fabs(z) == 0.0)
      { /* Map to center gray: */
        clr = (frgb_t){{ 0.500, 0.500, 0.500 }};
      }
    else
      { clr = frgb_path_map_signed_2(z, 1); }

    return clr;
  }

frgb_t fvw_texture_get_color
  ( fvw_texture_t *tx,
    int32_t x, int32_t y,
    float v_cur,
    float vmin_cur,
    float vmax_cur
  )
  {
    frgb_t clr;
    if (tx == NULL)
      { clr = fvw_texture_color_from_value(v_cur, vmin_cur, vmax_cur); }
    else
      { /* Get the texture image dimensions: */
        int32_t NC_tx, NX_tx, NY_tx;
        assert(tx != NULL);
        float_image_get_size(tx->timg, &NC_tx, &NX_tx, &NY_tx);
        assert((tx->chans[0] >= 0) && (tx->chans[0] < NC_tx));
        assert((x >= 0) && (x < NX_tx));
        assert((y >= 0) && (y < NY_tx));
        if (tx->chans[2] == -1)
          { assert(tx->chans[1] == -1);
            float v =  float_image_get_sample(tx->timg, tx->chans[0], x, y);
            clr = fvw_texture_color_from_value(v, tx->vmin, tx->vmax);
          }
        if (tx->chans[2] != -1)
          { assert((tx->chans[1] >= 0) && (tx->chans[1] < NC_tx));
            assert((tx->chans[2] >= 0) && (tx->chans[1] < NC_tx));
            for (int32_t k = 0; k < 3; k++)
              { clr.c[k] = float_image_get_sample(tx->timg, tx->chans[k], x, y); }
          }
      }
   return clr;
  }
