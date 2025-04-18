#ifndef fvw_paint_H
#define fvw_paint_H

/* fvw_paint.h - general painting tools. */
/* Last edited on 2025-04-14 05:12:33 by stolfi */

#include <stdint.h>

#define fvw_debug_GL FALSE
  /* If TRUE, report calls to the GL event methods. */

#define fvw_debug_paint FALSE
  /* If TRUE, report calls to internal painting routines. */

void fvw_compute_normal
  ( double xu, double yu, double zu,
    double xv, double yv, double zv,
    float *xnP, float *ynP, float *znP
  );
  /* Computes the normal {*xnP,*ynP,*znP} to a plane parallel to 
   the two vectors {(xu,yu,zu)} and {(xv,yv,zv)}. */
   
void fvw_paint_triangle
  ( float xa, float ya, float za,
    float xb, float yb, float zb,
    float xc, float yc, float zc
  );
  /* Paints the triangle whose corners are {(xa,ya,za),(xb,yb,zb),(xc,yc,zc)}. Assumes that
    the color has been set, and that the triangle is CCW as seen from above.
    Should be called between {glBegin(GL_TRIANGLES)} and {glEnd()}. */
   
void fvw_paint_quad
  ( float xa, float ya, float za,
    float xb, float yb, float zb,
    float xc, float yc, float zc,
    float xd, float yd, float zd
  );
  /* Paints the quadrilateral whose corners are {(xa,ya,za),(xb,yb,zb),(xc,yc,zc),(xd,yd,zd)}.
    Assumes that the color has been set.  Should be called between {glBegin(GL_QUADS)} and {glEnd()}. */

void fvw_paint_cell
  ( int32_t x, int32_t y, 
    float z00, float z10, float z01, float z11, 
    float CR, float CG, float CB
  );
  /* Paints the (non-planar) square patch of the terrain
    whose projection has opposite corners {(x,y)} and {(x+1,y+1)},
    with color {CR,CG,CB}. 
    
    The patch consists of four triangles with a common vertex ar
    {(xc,yc)=(x+0.5,y+0.5)} and the other vertices at {(x+dx,y+dy)}
    where {dx,dy} are 0 or 1. Assumes that the Z coordinates at the
    corners are {z00,z10,z01,z11}, respectively {SW,SE,NW,NE} (where +X
    is East and +Y is North). The height at {(xc,yc)} is assumed to be
    the mean of those four heights.
    
    Should be called between {glBegin(GL_TRIANGLES)} and {glEnd()}. */

void fvw_paint_node
  ( int32_t NX,
    int32_t NY,
    int32_t x,
    int32_t y,
    float zoo,
    float zom, float zpm, float zpo, float zpp, 
    float zop, float zmp, float zmo, float zmm,  
    float CR, float CG, float CB
  );
  /* Paints the (non-planar) unit square patch of the terrain centered
    on the corner {(x,y)} with color {CR,CG,CB}. 
    
    The patch consists of eight triangles whose vertices are
    {(x+dx,y+dy)}, for {dx,dy} in {-0.5,0,+0.5}, that share a common
    corner {(x,y)}. Assumes that the height at those nine points
    {(x+dx,y+dy)}, are {zoo,zom,...,zmm}; where the last two letters are
    the increments {dx} and {dy}, with {-0.5}, 0, and {+0.5} encoded as
    'm', 'o', 'p', respectively.
    
    Should be called between {glBegin(GL_TRIANGLES)} and {glEnd()}. */

void fvw_paint_horz_square
  ( float xlo, float xhi, 
    float ylo, float yhi, 
    float v,
    double zscale,
    float CR, float CG, float CB
  );
  /* Paints an horizontal square with X and Y ranges {[xlo_xhi]} and {[ylo_yhi]},
    with height {zscale*v} and color {CR,CG,CB}. 
    
    Should be called between {glBegin(GL_QUADS)} and {glEnd()}. */

void fvw_paint_vert_square
  ( float xa, float ya, 
    float xb, float yb, 
    float vlo, float vhi,
    double zscale,
    float CR, float CG, float CB
  );
  /* Paints a vertical square whose XY projection is the segment {(xa,ya)} to
    {(xb,yb)} with bottom Z at {zscale*vlo} and top Z at {zscale*vhi}
    and color {CR,CG,CB}.
    
    Should be called between {glBegin(GL_QUADS)} and {glEnd()}. */

#endif
