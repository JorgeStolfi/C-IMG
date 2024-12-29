#ifndef fvw_paint_cell_colored_H
#define fvw_paint_cell_colored_H

/* fvw_paint_cell_colored.h - painting a grid terrain with texmap-colored cells. */
/* Last edited on 2024-12-23 09:03:24 by stolfi */

#include <stdint.h>
#include <float_image.h>

void fvw_paint_cell_colored_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale,
    float_image_t *tx
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {zscale}.
    The texture {tx} must be non-null; paints each pixel between {(x,y)} and
    {(x+1,y+1)} with color {tx[x,y]} */

void fvw_paint_cell_colored_cell
  ( int32_t x, int32_t y, 
    float z00, float z10, float z01, float z11, 
    float CR, float CG, float CB
  );
  /* Paints the pixel whose lower left corner is {(x,y)}. Assumes that
    the Z coordinates at the corners are {z00,z10,z01,z11},
    respectively {SW,SE,NW,NE} (where +X is East and +Y is North).
    The pixel is painted with color {CR,CG,CB}, */



#endif 
