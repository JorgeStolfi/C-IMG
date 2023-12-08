/* See sample.h */
/* Last edited on 2021-08-20 17:20:16 by stolfi */

#include <math.h>
#include <affirm.h>

#include <imagebuf.h>
#include <sample.h>

double filter ARGS((double x));

void sample_pixel (ImageBuffer *im, double xp, double yp, int chan, double *pp)
{
  int dxmin, dxmax, dx, dymin, dymax, dy; 
  int y, x, yctr, xctr, yix;
  double xfr, yfr, xw, xwt, yw, ywt, xywt;
  double sw, swp;
  double *p;
  int nchans = im->nchans;
  
  /* fprintf(stderr, "sampling channel %d at (%.3f,%.3f)...\n", chan, xp,yp); */

  swp = 0.0; sw = 0.0;
  if (
      (yp >= (-(double)FILTER_RADIUS)) && (yp <= ((double)(im->rows-1 + FILTER_RADIUS))) &&
      (xp >= (-(double)FILTER_RADIUS)) && (xp <= ((double)(im->cols-1 + FILTER_RADIUS)))
    )
    { yctr = (int)(yp + FILTER_RADIUS) - FILTER_RADIUS; yfr = yp - (double)yctr;
      xctr = (int)(xp + FILTER_RADIUS) - FILTER_RADIUS; xfr = xp - (double)xctr;

      dymin = -FILTER_RADIUS; 
      if (yctr+dymin < 0) { dymin = -yctr; }
      dymax = FILTER_RADIUS;
      if (yctr+dymax > im->rows-1) { dymax = im->rows-1 - yctr; }

      for (dy = dymin; dy <= dymax; dy++)
        { y = yctr + dy;
          /* fprintf(stderr, " y = %d bufrows = %d", y, im->bufrows); */
          affirm((y >= im->inirow) & (y <= im->finrow), "bad row index");
          yix = ((y + FILTER_RADIUS*im->bufrows) % im->bufrows);
          yw = (double)dy - yfr;
          ywt = filter(yw);
          dxmin = -FILTER_RADIUS; 
          if (xctr+dxmin < 0) { dxmin = -xctr; }
          dxmax = FILTER_RADIUS; 
          if (xctr+dxmax > im->cols-1) { dxmax = im->cols-1 - xctr; }
          p = im->rowptr[yix];
          affirm(p != NULL, "oops - missing line");
          p += nchans*(xctr+dxmin) + chan;
          /* fprintf(stderr, " y = %d yix = %d p = %x", y, yix, (unsigned int)p); */
          for (dx = dxmin; dx <= dxmax; dx++)
            { x = xctr + dx;
              xw = (double)dx - xfr;
              xwt = filter(xw);
              xywt = xwt * ywt;
              /* fprintf(stderr, " [%d,%d]", y,x); */
              sw += xywt;
              swp += xywt * (*p);
              p += nchans;
            }
          /* fprintf(stderr, " p = %x\n", (unsigned int)p); */
        }
    }
  (*pp) = (sw != 0.0 ? swp/sw : 0.0);
  /* fprintf(stderr, "result = %.3f\n", *pp); */
}

double filter(double x)
{
  if (x == 0.0) 
    { return(1.0); }
  else
    { double z = (x/FILTER_SIGMA);
      return exp(-0.5*z*z)*sin(M_PI*x)/x;
    }
}

