#ifndef fvw_paint_H
#define fvw_paint_H

/* fvw_paint.h - general painting tools. */
/* Last edited on 2024-12-21 14:01:20 by stolfi */

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

#endif
