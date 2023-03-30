#ifndef fvw_paint_self_colored_hist_H
#define fvw_paint_self_colored_hist_H

/* fvw_paint_self_colored.h - painting a grid terrain without texture. */
/* Last edited on 2010-07-02 12:58:53 by stolfilocal */

#define _GNU_SOURCE
#include <float_image.h>

void fvw_paint_self_colored_hist_height_map
  ( float_image_t *ht, 
    int c, 
    double zscale, 
    float vmin,
    float vmax,
    bool_t skirt
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {zscale}.
    Each pixel is painted as a horizontal square, with vertical walls
    between adjacent pixels.
    The square's color is based on the sample value relative to the nominal sample 
    range {[vmin _ vmax]}.  If {skirt} is true paints the rod walls around the perimeter
    of the domain too, if false paints only the rod walls that are strictly inside the domain. */

void fvw_paint_self_colored_horz_square
  ( float xlo, float xhi, 
    float ylo, float yhi, 
    float v,
    double zscale,
    float vmin,
    float vmax 
  );
  /* Paints an horizontal square with X and Y ranges {[xlo_xhi]} and {[ylo_yhi]},
    with height {zscale*v}. The color is determined  by {v,vmin,vmax}. 
    Should be called between {glBegin(GL_QUADS)} and {glEnd()}. */

void fvw_paint_self_colored_vert_square
  ( float xa, float ya, 
    float xb, float yb, 
    float vlo, float vhi,
    double zscale,
    float vmin,
    float vmax 
  );
  /* Paints a vertical square whose XY projection is the segment {(xa,ya)} to
    {(xb,yb)} with bottom Z at {zscale*vlo} and top Z at {zscale*vhi}.
    The color is determined by {(vlo+vhi)/2,vmin,vmax}.
    Should be called between {glBegin(GL_QUADS)} and {glEnd()}. */

#endif
