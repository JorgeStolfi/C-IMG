#ifndef fvw_paint_height_map_H
#define fvw_paint_height_map_H

/* fvw_paint_height_map.h - painting a grid terrain with texmap-colored cells. */
/* Last edited on 2025-01-21 15:39:38 by stolfi */

#include <stdint.h>

#include <bool.h>
#include <float_image.h>
#include <fvw_texture.h>

void fvw_paint_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double ht_scale, 
    bool_t hist,
    fvw_texture_t *tx,
    float vmin_cur,
    float vmax_cur
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {ht_scale}.

    If {hist} is false, paints channel {c} of {ht} the image as a
    terrain, obtaining the color from {tx} and the other parameters as decribed
    under {fvw_texture_get_color}.
    
    If {hist} is true, paints channel {c} of {ht} as a histogram. In
    this case, if {tx.cell} is false, obtains the colors from {tx} and
    the other parameters as described under {fvw_texture_get_color}
    
    If {tx} is NULL, or {hist} is true and {tx.cell} is true, obtains
    the colors by converting each (unscaled) height value {v} itself to
    an RGB color with
    {fvw_texture_color_from_value(v,vmin_cur,vmax_cur)}. */

void fvw_paint_height_map_cell_colored
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale,
    fvw_texture_t *tx, 
    float vmin_cur, 
    float vmax_cur
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {zscale}.
    
    The flag {tx.cell} must be true, meaning that {tx.tmg} must be
    non-null and must have one col and one row less than {ht}. For each
    pixel {x,y} of {tx}, paints the unit square of the terrain with
    corners {(x,y)} and {(x+1,y+1)} with color defined by {tx} and the
    other parameters, as per {fvw_texture_get_color}. */

void fvw_paint_height_map_node_colored
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale,
    fvw_texture_t *tx, 
    float vmin_cur, 
    float vmax_cur
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {zscale}.
    
    The flag {tx.cell} must be false, meaning that the texture image
    {tx.timg} has the same size as {ht}. For each pixel {x,y} of {ht},
    paints the unit square of the terrain centered on the corner {(x,y)}
    with color defined by {tx} and the other parameters, as per
    {fvw_texture_get_color}. */

void fvw_paint_height_map_hist_colored
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale,
    fvw_texture_t *tx, 
    float vmin_cur,
    float vmax_cur,
    bool_t skirt
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {zscale}.
    Each pixel is painted as a horizontal square, with vertical walls
    between adjacent pixels.
    
    If {tx} is not {NULL}, then the flag {tx.cell} must be false.
    The procedure then paints the unit square of the
    terrain centered on each corner {(x,y)} with color defined by {tx}
    and the other parameters, as per {fvw_texture_get_color}. 
    
    If {tx} is {NULL}, and obtains the colors by converting
    each (unscaled) height value {v} to an RGB color with
    {fvw_texture_color_from_value(v,vmin_cur,vmax_cur)}.
    
    If {skirt} is true paints the rod walls around the perimeter of the
    domain too, if false paints only the rod walls that are strictly
    inside the domain. */
 
#endif 
