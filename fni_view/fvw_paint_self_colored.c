/* See fvw_paint_self_colored.h */
/* Last edited on 2024-12-23 09:06:33 by stolfi */

#include <assert.h>
#include <math.h>
#include <GL/glu.h>

#include <float_image.h>
#include <frgb.h>
#include <frgb_path.h>

#include <fvw_paint.h>
#include <fvw_paint_self_colored.h>

void fvw_paint_self_colored_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale, 
    float vmin,
    float vmax
  )
  {
    if (fvw_debug_paint) { fprintf(stderr, "+ %s\n", __FUNCTION__); }
    
    /* Get  the height image dimensons: */
    assert(ht != NULL);
    int32_t HNC, HNX, HNY;
    float_image_get_size(ht, &HNC, &HNX, &HNY);

    /* We need buffers for two rows of height map samples: */
    float va[HNX], vb[HNX];
    float *v0 = va; /* Current row of height samples. */
    float *v1 = vb; /* Next row of height samples. */
    
    /* We need buffers for two rows of mid-vert-edge height means: */
    float ma[HNX], mb[HNX];
    float *m0 = ma; /* Height sample means for row {y-1/2}. */
    float *m1 = mb; /* Height sample means for row {y+1/2}. */
    
    /* Get first row of samples: */
    float_image_get_sample_row(ht, (int32_t)c, 0, HNX-1, 0, v0);
    
    /* Scan rows of height array: */
    int32_t x, y;
    for(y = 0; y < HNY; y++)
      { if (y < HNY-1)
          { /* Get in {v1} the heights at ordinate {y+1}: */
            float_image_get_sample_row(ht, (int32_t)c, 0, HNX-1, y+1, v1);
            /* Compute the means {m1} for row {y+1/2}: */
            for(x = 0; x < HNX; x++)
              { /* Compute the mean height for vertical edge at {x}: */
                m1[x] = (v0[x] + v1[x])/2;
              }
          }
        /* Now paint 8-triangle patches around corners in row {y}: */
        float *pm0 = m0; /* Pointer to mid-edge height in col {x}, row {y-1/2}. */
        float *pm1 = m1; /* Pointer to mid-edge height in col {x}, row {y+1/2}. */
        float *pv0 = v0; /* Pointer to corner height in col {x}, row {y}. */
        for(x = 0; x < HNX; x++)
          { glBegin(GL_TRIANGLES);
            /* Get the height {voo} at {(x,y)}: */ 
            float xo = (float)x, xm = (float)(x - 0.5), xp = (float)(x + 0.5);
            float yo = (float)y, ym = (float)(y - 0.5), yp = (float)(y + 0.5);
            float voo = pv0[0];
            /* Get the mid-edge heights around {(x,y)}: */
            float vmo = (x > 0 ? (pv0[-1]+pv0[0])/2 : 0);
            float vpo = (x < HNX-1 ? (pv0[0]+pv0[1])/2 : 0);
            float vom = (y > 0 ? pm0[0] : 0);
            float vop = (y < HNY-1 ? pm1[0] : 0);
            if (y > 0)
              { /* Paint the four triangles between ordinates {y-1/2} and {y}: */
                if (x > 0)
                  { float vmm = (pm0[-1] + pm0[0])/2;
                    fvw_paint_self_colored_triangle(xo, yo, voo, xm, yo, vmo, xm, ym, vmm, zscale, vmin, vmax);
                    fvw_paint_self_colored_triangle(xo, yo, voo, xm, ym, vmm, xo, ym, vom, zscale, vmin, vmax);
                  }
                if (x < HNX-1)
                  { float vpm = (pm0[0] + pm0[1])/2;
                    fvw_paint_self_colored_triangle(xo, yo, voo, xo, ym, vom, xp, ym, vpm, zscale, vmin, vmax);
                    fvw_paint_self_colored_triangle(xo, yo, voo, xp, ym, vpm, xp, yo, vpo, zscale, vmin, vmax);
                  }
              }
            if (y < HNY-1)
              { /* Paint the four triangles between ordinates {y} and {y+1/2}: */
                if (x < HNX-1)
                  { float vpp = (pm1[0] + pm1[1])/2;
                    fvw_paint_self_colored_triangle(xo, yo, voo, xp, yo, vpo, xp, yp, vpp, zscale, vmin, vmax);
                    fvw_paint_self_colored_triangle(xo, yo, voo, xp, yp, vpp, xo, yp, vop, zscale, vmin, vmax);
                  }
                if (x > 0)
                  { float vmp = (pm1[-1] + pm1[0])/2;
                    fvw_paint_self_colored_triangle(xo, yo, voo, xo, yp, vop, xm, yp, vmp, zscale, vmin, vmax);
                    fvw_paint_self_colored_triangle(xo, yo, voo, xm, yp, vmp, xm, yo, vmo, zscale, vmin, vmax);
                  }
              }
            glEnd();
            pv0++;
            pm0++;
            pm1++;
          }
        /* Swap row buffers, {v1-->v0}, {m1<-->m0}: */
        { float *t = v0; v0 = v1; v1 = t; }
        { float *t = m0; m0 = m1; m1 = t; }
      }
    if (fvw_debug_paint) { fprintf(stderr, "- %s\n", __FUNCTION__); }
  }
  
void fvw_paint_self_colored_triangle
  ( float xa, float ya, float va,
    float xb, float yb, float vb,
    float xc, float yc, float vc,
    double zscale,
    float vmin,
    float vmax 
  )
  {
    /* Compute height at barycenter of triangle: */
    float vm = (float)((va + vb + vc)/3.0);
    
    /* Get {vzer,vdel} so that {vzer+vdel}->{cmax}, {vzer-vdel}->complement. */ 
    double vzer, vdel;
    double dv = vmax - vmin;
    double eps = 1.0e-4*dv;
    if (vmin > 0)
      { vzer = vmin - eps; vdel = dv + 2*eps; }
    else if (vmax < 0)
      { vzer = vmax + eps; vdel = dv + 2*eps; }
    else 
      { vzer = 0.0; vdel = fmax(vmax, -vmin); }

    /* Make sure that {vdel} is nonzero: */
    if (vdel == 0) { vdel = 1.0; }
    
    float clr[3];
    fvw_color_from_value(3, vm - vzer, vdel, clr);
    GLfloat CR = clr[0], CG = clr[1], CB = clr[2];
    glColor3f(CR,CG,CB);
    float sva = (float)(zscale*va);
    float svb = (float)(zscale*vb);
    float svc = (float)(zscale*vc);
    fvw_paint_triangle(xa, ya, sva, xb, yb, svb, xc, yc, svc);
  }

#define Ymax (0.8667)
  /* Max brightness to use when painting in grays. */

#define Ymin (0.5333)
  /* Min brightness for positive values, when painting in grays. */

void fvw_color_from_value(int32_t NC, double v, double vdel, float clr[])
  {
    /* Compute perceptual brightness {z} in {[0_1]}: */
    double z;
    if (v < -vdel)
      { z = -1.0; }
    else if (v > +vdel)
      { z = +1.0; }
    else if (vdel == 0)
      { z = 0.5; }
    else 
      { z = v/vdel; }
      
    uint32_t c;
    if (fabs(z) == 0.0)
      { /* Map to center gray: */
        for (c = 0; c < NC; c++) { clr[c] = 0.500; }
      }
    else if (NC == 1)
      { /* The luminance interpolates from {ymax} to {ymin} with {abs(z)} as ratio: */
        double az = fabs(z);
        double ymin = 0.550;
        double ymax = 0.900;
        double smax = az, smin = 1.0 - az;
        double y = smin*ymin + smax*ymax; 
        /* If {z} is negative, complement the color relative to middle gray: */
        if (z < 0) { y = 1.000 - y; }
        clr[0] = (float)y;
      }
    else
      { frgb_t fv = frgb_path_map_signed_2(z, 1);
        for (c = 0; c < NC; c++) { clr[c] = fv.c[c]; }
      }
  }

