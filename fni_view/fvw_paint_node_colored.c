/* See fvw_paint_node_colored.h */
/* Last edited on 2017-06-25 16:26:35 by stolfilocal */

#define _GNU_SOURCE
#include <assert.h>
#include <GL/glu.h>

#include <float_image.h>

#include <fvw_paint.h>
#include <fvw_paint_node_colored.h>

void fvw_paint_node_colored_height_map
  ( float_image_t *ht, 
    int c, 
    double zscale,
    float_image_t *tx
  )
  {
    /* Get  the height image dimensons: */
    assert(ht != NULL);
    int HNC, HNX, HNY;
    float_image_get_size(ht, &HNC, &HNX, &HNY);
    
    /* Get the texture image dimensions: */
    int TNC, TNX, TNY;
    assert(tx != NULL);
    float_image_get_size(tx, &TNC, &TNX, &TNY);
    assert((TNC == 1) || (TNC == 3));
    assert(TNX == HNX);
    assert(TNY == HNY);

    /* We need buffers for two rows of height map samples: */
    float va[HNX], vb[HNX];
    float *v0 = va; /* Current row of height samples. */
    float *v1 = vb; /* Next row of height samples. */
    
    /* We need buffers for two rows of mid-vert-edge height means: */
    float ma[HNX], mb[HNX];
    float *m0 = ma; /* Height sample means for row {y-1/2}. */
    float *m1 = mb; /* Height sample means for row {y+1/2}. */

    /* We need a buffer for one row of texture colors: */
    float clr[TNC*TNX];

    /* Get in {v0} the heights of row 0: */
    float_image_get_sample_row(ht, c, 0, HNX-1, 0, v0);
    
    /* Scan rows of height array: */
    int x, y;
    for(y = 0; y < HNY; y++)
      { /* Get in {clr} the colors of corners at ordinate {y}: */
        float_image_get_pixel_row(tx, 0, TNX-1, y, clr);
        if (y < HNY-1)
          { /* Get in {v1} the heights at ordinate {y+1}: */
            float_image_get_sample_row(ht, c, 0, HNX-1, y+1, v1);
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
        float *pc = clr; /* Pointer to color or corner in col {x}, row {y}. */
        for(x = 0; x < HNX; x++)
          { /* Get color in GL format: */
            GLfloat CR, CG, CB;
            if (TNC == 1)
              { CR = CG = CB = pc[0]; }
            else
              { CR = pc[0]; CG = pc[1]; CB = pc[2]; }
            pc += TNC;
            glColor3f(CR,CG,CB);
            /* !!! Should use {GL_TRIANGLE_FAN} instead of {GL_TRIANGLES}. !!! */
            glBegin(GL_TRIANGLES);
            /* Get the height {zoo} at {(x,y)}: */ 
            float xo = (float)x, xm = (float)(x - 0.5), xp = (float)(x + 0.5);
            float yo = (float)y, ym = (float)(y - 0.5), yp = (float)(y + 0.5);
            float zoo = (float)(zscale*pv0[0]);
            /* Get the mid-edge heights around {(x,y)}: */
            float zmo = (float)(x > 0 ? zscale*(pv0[-1]+pv0[0])/2 : 0);
            float zpo = (float)(x < HNX-1 ? zscale*(pv0[0]+pv0[1])/2 : 0);
            float zom = (float)(y > 0 ? zscale*pm0[0] : 0);
            float zop = (float)(y < HNY-1 ? zscale*pm1[0] : 0);
            if (y > 0)
              { /* Paint the four triangles between ordinates {y-1/2} and {y}: */
                if (x > 0)
                  { float zmm = (float)(zscale*(pm0[-1] + pm0[0])/2);
                    fvw_paint_triangle(xo, yo, zoo, xm, yo, zmo, xm, ym, zmm);
                    fvw_paint_triangle(xo, yo, zoo, xm, ym, zmm, xo, ym, zom);
                  }
                if (x < HNX-1)
                  { float zpm = (float)(zscale*(pm0[0] + pm0[1])/2);
                    fvw_paint_triangle(xo, yo, zoo, xo, ym, zom, xp, ym, zpm);
                    fvw_paint_triangle(xo, yo, zoo, xp, ym, zpm, xp, yo, zpo);
                  }
              }
            if (y < HNY-1)
              { /* Paint the four triangles between ordinates {y} and {y+1/2}: */
                if (x < HNX-1)
                  { float zpp = (float)(zscale*(pm1[0] + pm1[1])/2);
                    fvw_paint_triangle(xo, yo, zoo, xp, yo, zpo, xp, yp, zpp);
                    fvw_paint_triangle(xo, yo, zoo, xp, yp, zpp, xo, yp, zop);
                  }
                if (x > 0)
                  { float zmp = (float)(zscale*(pm1[-1] + pm1[0])/2);
                    fvw_paint_triangle(xo, yo, zoo, xo, yp, zop, xm, yp, zmp);
                    fvw_paint_triangle(xo, yo, zoo, xm, yp, zmp, xm, yo, zmo);
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
  }
