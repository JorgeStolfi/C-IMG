#ifndef fvw_paint_node_colored_H
#define fvw_paint_node_colored_H

/* fvw_paint_node_colored.h - painting a grid terrain with texmap-colored vertices. */
/* Last edited on 2024-12-23 09:02:32 by stolfi */

#include <stdint.h>
#include <float_image.h>

void fvw_paint_node_colored_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale,
    float_image_t *tx
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {zscale}.
    The image {tx} must be non-null; paints each pixel between {(x,y)} and
    {(x+1,y+1)} with color {tx[x,y]}. */
 
#endif
