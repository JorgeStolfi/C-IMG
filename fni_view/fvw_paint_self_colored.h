#ifndef fvw_paint_self_colored_H
#define fvw_paint_self_colored_H

/* fvw_paint_self_colored.h - painting a grid terrain without texture. */
/* Last edited on 2010-07-02 12:58:53 by stolfilocal */

#define _GNU_SOURCE
#include <float_image.h>

void fvw_paint_self_colored_height_map
  ( float_image_t *ht, 
    int c, 
    double zscale, 
    float vmin,
    float vmax
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {zscale}.
    Each pixel is divided into 4 triangles, and each is painted
    with a color based on the sample value relative to the nominal sample 
    range {[vmin _ vmax]}. */

void fvw_paint_self_colored_triangle
  ( float xa, float ya, float va,
    float xb, float yb, float vb,
    float xc, float yc, float vc,
    double zscale,
    float vmin,
    float vmax 
  );
  /* Paints the triangle with corners {(xa,ya,zscale*va)}, {(xb,yb,zscale*vb)}, {(xc,yc,zscale*vc)}.
    The color depends on the position of {va,vb,vc} relative to the range {[vmin_vmax]}.
    Assumes that the triangle is CCW as seen from above.
    Should be called between {glBegin(GL_TRIANGLES)} and {glEnd()}. */

void fvw_color_from_value(int NC, double v, double vdel, float clr[]);
  /* Computes artificial color {clr[0..NC-1]} for a pixel of value {v},
     whose nominal range is {[-vdel _ +vdel]}. */

#endif
