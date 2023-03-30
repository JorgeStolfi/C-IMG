/* See fvw_paint.h */
/* Last edited on 2013-10-21 03:33:39 by stolfilocal */

#define _GNU_SOURCE
#include <math.h>
#include <GL/glu.h>

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
