/* See fvw_paint_cell_colored.h */
/* Last edited on 2024-12-23 09:06:15 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <GL/glu.h>

#include <float_image.h>

#include <fvw_paint.h>
#include <fvw_paint_cell_colored.h>

void fvw_paint_cell_colored_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale,
    float_image_t *tx
  )
  {
    /* Get  the height image dimensons: */
    assert(ht != NULL);
    int32_t HNC, HNX, HNY;
    float_image_get_size(ht, &HNC, &HNX, &HNY);
    
     /* Get the texture image dimensions: */
    int32_t TNC, TNX, TNY;
    assert(tx != NULL);
    float_image_get_size(tx, &TNC, &TNX, &TNY);
    assert((TNC == 1) || (TNC == 3));
    assert(TNX == HNX-1);
    assert(TNY == HNY-1);

    /* We need buffers for two rows of height map samples: */
    float va[HNX], vb[HNX];
    float *v0 = va; /* Current row of height samples. */
    float *v1 = vb; /* Next row of height samples. */
    
    /* We need a buffer for one row of texture colors: */
    float clr[TNC*TNX];

    /* Get first row of samples: */
    float_image_get_sample_row(ht, (int32_t)c, 0, HNX-1, 0, v0);
    
    /* Scan rows of height array: */
    int32_t x, y;
    for(y = 0; y < HNY-1; y++)
      { /* Get next row of samples: */
        float_image_get_sample_row(ht, (int32_t)c, 0, HNX-1, y+1, v1);
        /* Get next row of colors: */
        float_image_get_pixel_row(tx, 0, TNX-1, y, clr);
        /* Now paint pixels: */
        float *pv0 = v0; /* Pointer to height in col {x}, row {y}. */
        float *pv1 = v1; /* Pointer to height in col {x}, row {y+1}. */
        float *pc = clr; /* Pointer to color of pixel in col {x}, row {y}. */
        for(x = 0; x < HNX-1; x++)
          { /* Get color in GL format: */
            GLfloat CR, CG, CB;
            if (TNC == 1)
              { CR = CG = CB = pc[0]; }
            else
              { CR = pc[0]; CG = pc[1]; CB = pc[2]; }
            pc += TNC;
            /* Get heights at four pixel corners in GL format: */
            float z00 = (float)(zscale*pv0[0]);  /* Height at {(x,   y  )}. */
            float z10 = (float)(zscale*pv0[1]);  /* Height at {(x+1, y  )}. */
            float z01 = (float)(zscale*pv1[0]);  /* Height at {(x,   y+1)}. */
            float z11 = (float)(zscale*pv1[1]);  /* Height at {(x+1, y+1)}. */
            fvw_paint_cell_colored_cell(x,y, z00,z10,z01,z11, CR,CG,CB);
            pv0++; pv1++;
          }
        /* Swap row buffers, {v1} now becomes {v0}: */
        { float *t = v0; v0 = v1; v1 = t; }
      }
  }
            
void fvw_paint_cell_colored_cell
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
