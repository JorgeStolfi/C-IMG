/* See fvw_paint.h */
/* Last edited on 2025-04-14 05:12:34 by stolfi */

#include <stdint.h>
#include <math.h>
#include <GL/glu.h>

#include <fvw_texture.h>

#include <fvw_paint.h>

void fvw_compute_normal
  ( double xu, double yu, double zu,
    double xv, double yv, double zv,
    float *xnP, float *ynP, float *znP
  )
  {
    double xn = yu*zv - yv*zu;
    double yn = zu*xv - zv*xu;
    double zn = xu*yv - xv*yu;
    double nmod = sqrt(xn*xn + yn*yn + zn*zn);
    (*xnP)= (float)(xn / nmod);
    (*ynP)= (float)(yn / nmod);
    (*znP)= (float)(zn / nmod);
  }

void fvw_paint_triangle
  ( float xa, float ya, float za,
    float xb, float yb, float zb,
    float xc, float yc, float zc
  )
  {
    double xab = xb-(double)xa;
    double yab = yb-(double)ya;
    double zab = zb-(double)za;
    
    double xac = xc-(double)xa;
    double yac = yc-(double)ya;
    double zac = zc-(double)za;
    
    float xn, yn, zn;
    fvw_compute_normal(xab, yab, zab, xac, yac, zac, &xn, &yn, &zn);

    glNormal3f(xn, yn, zn);

    glVertex3f(xa, ya, za);
    glVertex3f(xb, yb, zb);
    glVertex3f(xc, yc, zc);
  }

void fvw_paint_quad
  ( float xa, float ya, float za,
    float xb, float yb, float zb,
    float xc, float yc, float zc,
    float xd, float yd, float zd
  )
  {
    /* Compute normal assuming that it is flat: */
    double xab = xb-(double)xa;
    double yab = yb-(double)ya;
    double zab = zb-(double)za;
    
    double xad = xd-(double)xa;
    double yad = yd-(double)ya;
    double zad = zd-(double)za;
    
    float xn, yn, zn;
    fvw_compute_normal(xab, yab, zab, xad, yad, zad, &xn, &yn, &zn);

    glNormal3f(xn, yn, zn);

    glVertex3f(xa, ya, za);
    glVertex3f(xb, yb, zb);
    glVertex3f(xc, yc, zc);
    glVertex3f(xd, yd, zd);
  }
            
void fvw_paint_cell
  ( int32_t x, int32_t y, 
    float z00, float z10, float z01, float z11, 
    float CR, float CG, float CB
  )
  {
    /* Low and high coordinates of pixel on XY plane: */
    float x0 = (float)x;
    float y0 = (float)y;
    float x1 = x0 + 1.00f;
    float y1 = y0 + 1.00f;
    
    glColor3f(CR,CG,CB);
    /*
      Plot pixel as four triangles {p,q,u} where {u} is the center 
      and {p,q} scan the boundary in the order
        01<----11
        |       ^
        |  .u   |
        V       |
        00---->10
    */

    /* The height at the center is the average of the four corner heights: */
    float xu = (float)(x + 0.5);
    float yu = (float)(y + 0.5);
    float zu = (z00 + z10 + z01 + z11)/4;

    /* !!! Should use {GL_TRIANGLE_FAN} instead of {GL_TRIANGLES}. !!! */
    glBegin(GL_TRIANGLES);
    fvw_paint_triangle(xu, yu, zu, x0, y0, z00, x1, y0, z10);
    fvw_paint_triangle(xu, yu, zu, x1, y0, z10, x1, y1, z11);
    fvw_paint_triangle(xu, yu, zu, x1, y1, z11, x0, y1, z01);
    fvw_paint_triangle(xu, yu, zu, x0, y1, z01, x0, y0, z00);
    glEnd();
  }

void fvw_paint_node
  ( int32_t NX,
    int32_t NY,
    int32_t x,
    int32_t y,
    float zoo,
    float zom, float zpm, float zpo,
    float zpp, float zop, float zmp,
    float zmo, float zmm,  
    float CR, float CG, float CB
  )
  {
    float xo = (float)x, xm = (float)(x - 0.5), xp = (float)(x + 0.5);
    float yo = (float)y, ym = (float)(y - 0.5), yp = (float)(y + 0.5);
    if (y > 0)
      { /* Paint the four triangles between ordinates {y-1/2} and {y}: */
        if (x > 0)
          { fvw_paint_triangle(xo, yo, zoo, xm, yo, zmo, xm, ym, zmm);
            fvw_paint_triangle(xo, yo, zoo, xm, ym, zmm, xo, ym, zom);
          }
        if (x < NX-1)
          { fvw_paint_triangle(xo, yo, zoo, xo, ym, zom, xp, ym, zpm);
            fvw_paint_triangle(xo, yo, zoo, xp, ym, zpm, xp, yo, zpo);
          }
      }
    if (y < NY-1)
      { /* Paint the four triangles between ordinates {y} and {y+1/2}: */
        if (x < NX-1)
          { fvw_paint_triangle(xo, yo, zoo, xp, yo, zpo, xp, yp, zpp);
            fvw_paint_triangle(xo, yo, zoo, xp, yp, zpp, xo, yp, zop);
          }
        if (x > 0)
          { fvw_paint_triangle(xo, yo, zoo, xo, yp, zop, xm, yp, zmp);
            fvw_paint_triangle(xo, yo, zoo, xm, yp, zmp, xm, yo, zmo);
          }
      }
  }

void fvw_paint_horz_square
  ( float xlo, float xhi, 
    float ylo, float yhi, 
    float v,
    double zscale,
    float CR, float CG, float CB 
  )
  {
    glColor3f(CR,CG,CB);
    
    GLfloat z = (GLfloat)(zscale*v);
    fvw_paint_quad(xlo, ylo, z, xhi, ylo, z, xhi, yhi, z, xlo, yhi, z);
  }

void fvw_paint_vert_square
  ( float xa, float ya, 
    float xb, float yb, 
    float vlo, float vhi,
    double zscale,
    float CR, float CG, float CB
  )
  {
    glColor3f(CR,CG,CB);
    
    GLfloat zlo = (GLfloat)(zscale*vlo);
    GLfloat zhi = (GLfloat)(zscale*vhi);
    fvw_paint_quad(xa, ya, zlo, xb, yb, zlo, xb, yb, zhi, xa, ya, zhi);
  }
