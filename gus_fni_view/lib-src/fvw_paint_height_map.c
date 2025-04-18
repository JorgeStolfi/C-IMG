/* See fvw_paint_height_map.h */
/* Last edited on 2025-01-22 07:19:13 by stolfi */

#include <stdint.h>
#include <assert.h>
#include <GL/glu.h>

#include <bool.h>
#include <float_image.h>

#include <fvw_paint.h>
#include <fvw_texture.h>

#include <fvw_paint_height_map.h>

void fvw_paint_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double ht_scale, 
    bool_t hist,
    fvw_texture_t *tx,
    float vmin_cur,
    float vmax_cur
  )
  {
    if (ht == NULL) { return; }
    if (fvw_debug_paint) { fprintf(stderr, "+ %s\n", __FUNCTION__); }

    /* Get  the height image dimensons: */
    int32_t ht_NC, ht_NX, ht_NY;
    float_image_get_size(ht, &ht_NC, &ht_NX, &ht_NY);

    /* Set surface finish: */
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    
    if (hist)
      { /* Paint as histogram, no skirt: */
        fvw_paint_height_map_hist_colored(ht, c, ht_scale, tx, vmin_cur, vmax_cur, FALSE);
      }
    else if ((tx == NULL) || (! tx->cell))
      { /* TPaint as terrain, texture colors are associated with grid corners: */ 
        fvw_paint_height_map_node_colored(ht, c, ht_scale, tx, vmin_cur, vmax_cur);
      }
    else if ((tx != NULL) && tx->cell)
      { /* Paint as terrain, texture colors are associated with grid cells: */
        fvw_paint_height_map_cell_colored(ht, c, ht_scale, tx, vmin_cur, vmax_cur);
      }
    else
      { assert(FALSE); }
    glDisable(GL_COLOR_MATERIAL);
    if (fvw_debug_paint) { fprintf(stderr, "- %s\n", __FUNCTION__); }
  }

void fvw_paint_height_map_cell_colored
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale,
    fvw_texture_t *tx, 
    float vmin_cur, 
    float vmax_cur
  )
  {
    /* Get  the height image dimensons: */
    assert(ht != NULL);
    int32_t NC_ht, NX_ht, NY_ht;
    float_image_get_size(ht, &NC_ht, &NX_ht, &NY_ht);

    /* We need buffers for two rows of height map samples: */
    float va[NX_ht], vb[NX_ht];
    float *v0 = va; /* Current row of height samples. */
    float *v1 = vb; /* Next row of height samples. */
    
    /* Get first row of samples: */
    float_image_get_sample_row(ht, (int32_t)c, 0, NX_ht-1, 0, v0);
    
    /* Scan rows of height array: */
    int32_t x, y;
    for(y = 0; y < NY_ht-1; y++)
      { /* Get next row of samples: */
        float_image_get_sample_row(ht, (int32_t)c, 0, NX_ht-1, y+1, v1);
        /* Now paint pixels: */
        float *pv0 = v0; /* Pointer to height in col {x}, row {y}. */
        float *pv1 = v1; /* Pointer to height in col {x}, row {y+1}. */
        for(x = 0; x < NX_ht-1; x++)
          { /* Get unscaled heights at four pixel corners: */
            float v00 = pv0[0];  /* Unscaled height at {(x,   y  )}. */
            float v10 = pv0[1];  /* Unscaled height at {(x+1, y  )}. */
            float v01 = pv1[0];  /* Unscaled height at {(x,   y+1)}. */
            float v11 = pv1[1];  /* Unscaled height at {(x+1, y+1)}. */
            /* Get scaled heights at four pixel corners in GL format: */
            float z00 = (float)(zscale*v00);  /* Scaled height at {(x,   y  )}. */
            float z10 = (float)(zscale*v10);  /* Scaled height at {(x+1, y  )}. */
            float z01 = (float)(zscale*v01);  /* Scaled height at {(x,   y+1)}. */
            float z11 = (float)(zscale*v11);  /* Scaled height at {(x+1, y+1)}. */
            /* Get cell color: */
            frgb_t clr = fvw_texture_get_color(tx, x, y, v00, vmin_cur, vmax_cur);
            GLfloat CR = clr.c[0], CG = clr.c[1], CB = clr.c[2];
            fvw_paint_cell(x,y, z00,z10,z01,z11, CR,CG,CB);
            pv0++; pv1++;
          }
        /* Swap row buffers, {v1} now becomes {v0}: */
        { float *t = v0; v0 = v1; v1 = t; }
      }
  }

void fvw_paint_height_map_node_colored
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale,
    fvw_texture_t *tx, 
    float vmin_cur, 
    float vmax_cur
  )
  {
    /* Get  the height image dimensons: */
    assert(ht != NULL);
    int32_t NC_ht, NX_ht, NY_ht;
    float_image_get_size(ht, &NC_ht, &NX_ht, &NY_ht);

    /* We need buffers for two rows of height map samples: */
    float va[NX_ht], vb[NX_ht];
    float *v0 = va; /* Current row of height samples. */
    float *v1 = vb; /* Next row of height samples. */
    
    /* We need buffers for two rows of mid-vert-edge height means: */
    float ma[NX_ht], mb[NX_ht];
    float *m0 = ma; /* Height sample means for row {y-1/2}. */
    float *m1 = mb; /* Height sample means for row {y+1/2}. */

    /* Get in {v0} the heights of row 0: */
    float_image_get_sample_row(ht, (int32_t)c, 0, NX_ht-1, 0, v0);
    
    /* Scan rows of height array: */
    int32_t x, y;
    for(y = 0; y < NY_ht; y++)
      { if (y < NY_ht-1)
          { /* Get in {v1} the heights at ordinate {y+1}: */
            float_image_get_sample_row(ht, (int32_t)c, 0, NX_ht-1, y+1, v1);
            /* Compute the means {m1} for row {y+1/2}: */
            for(x = 0; x < NX_ht; x++)
              { /* Compute the mean height for vertical edge at {x}: */
                m1[x] = (v0[x] + v1[x])/2;
              }
          }
        /* Now paint 8-triangle patches around corners in row {y}: */
        float *pm0 = m0; /* Pointer to mid-edge height in col {x}, row {y-1/2}. */
        float *pm1 = m1; /* Pointer to mid-edge height in col {x}, row {y+1/2}. */
        float *pv0 = v0; /* Pointer to corner height in col {x}, row {y}. */
        for(x = 0; x < NX_ht; x++)
          { /* Get the unscaled height {voo} at {(x,y)}: */ 
            float voo = pv0[0];
            /* Get the scaled heights at {(x+dx,y+dy)}: */ 
            float zoo = (float)(zscale*voo);
            float zmo = (float)(x > 0 ? zscale*(pv0[-1]+pv0[0])/2 : 0);
            float zpo = (float)(x < NX_ht-1 ? zscale*(pv0[0]+pv0[1])/2 : 0);
            float zom = (float)(y > 0 ? zscale*pm0[0] : 0);
            float zop = (float)(y < NY_ht-1 ? zscale*pm1[0] : 0);
            float zmm = (float)((y > 0) && (x > 0) ? zscale*(pm0[-1] + pm0[0])/2 : 0);
            float zpm = (float)((y > 0) && (x < NX_ht-1) ? zscale*(pm0[0] + pm0[1])/2 : 0);
            float zpp = (float)((y < NY_ht-1) && (x < NX_ht-1) ? zscale*(pm1[0] + pm1[1])/2 : 0);
            float zmp = (float)((y < NY_ht-1) && (x > 0) ? zscale*(pm1[-1] + pm1[0])/2 : 0);
                        
            /* Get vertex color: */
            frgb_t clr = fvw_texture_get_color(tx, x, y, voo, vmin_cur, vmax_cur);
            GLfloat CR = clr.c[0], CG = clr.c[1], CB = clr.c[2];
            glColor3f(CR,CG,CB);
            /* !!! Should use {GL_TRIANGLE_FAN} instead of {GL_TRIANGLES}. !!! */
            glBegin(GL_TRIANGLES);
            fvw_paint_node(NX_ht, NY_ht, x,y, zoo, zom,zpm,zpo,zpp,zop,zmp,zmo,zmm, CR,CG,CB);
  
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

void fvw_paint_height_map_hist_colored
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale, 
    fvw_texture_t *tx, 
    float vmin_cur,
    float vmax_cur,
    bool_t skirt
  )
  {
    if (fvw_debug_paint) { fprintf(stderr, "+ %s\n", __FUNCTION__); }
    
    /* Get  the height image dimensons: */
    assert(ht != NULL);
    int32_t NC_ht, NX_ht, NY_ht;
    float_image_get_size(ht, &NC_ht, &NX_ht, &NY_ht);

    /* We need buffers for two rows of height map samples: */
    float va[NX_ht], vb[NX_ht];
    float *vo = va; /* Current row of height samples. */
    float *vm = vb; /* Next row of height samples. */
    
    /* Scan rows of height array: */
    int32_t x, y;
    for(y = 0; y <= NY_ht; y++)
      { /* If {y > 0}, then {vm} contains the pixel values of row {y-1}. */
	if (y <= NY_ht-1)
          { /* Get in {vo} the pixel values of row {y}: */
            float_image_get_sample_row(ht, (int32_t)c, 0, NX_ht-1, y, vo);
          }
	/* Paint the pixels in row {y} and the walls between row {y} and row {y-1}: */
	glBegin(GL_QUADS);
	float *pvo = vo; /* Pointer to pixel value in col {x}, row {y}. */
	float vmo = 0;   /* Pixel value in col {x-1}, row {y} */
        float *pvm = vm; /* Pointer to pixel value in col {x}, row {y-1}. */
	for(x = 0; x <= NX_ht; x++)
	   { float xlo = (float)x, xhi = (float)(x + 1.0); /* X range of rod. */
             float ylo = (float)y, yhi = (float)(y + 1.0); /* Y range of rof */
             float voo = (float)(x < NX_ht ? (*pvo) : 0); /* Pixel value in col {x}, row {y} */
             float vom = (float)(x < NX_ht ? (*pvm) : 0); /* Pixel value in col {x}, row {y-1} */
             /* Paint top of rod: */
	     if(( x < NX_ht) && ( y < NY_ht))
	       { frgb_t clr = fvw_texture_get_color(tx, x, y, voo, vmin_cur, vmax_cur);
                 fvw_paint_horz_square(xlo, xhi, ylo, yhi, voo, zscale, clr.c[0], clr.c[1], clr.c[2]); 
               }
             /* Paint west wall of rod: */
             if ((y < NY_ht) && (skirt || ((x > 0) && (x < NX_ht))))
               { frgb_t clr = fvw_texture_get_color(tx, x, y, (vmo+voo)/2, vmin_cur, vmax_cur);
                 fvw_paint_vert_square(xlo, ylo, xlo, yhi, vmo, voo, zscale, clr.c[0], clr.c[1], clr.c[2]);
               }
             /* Paint south wall of rod: */
             if (skirt || ((y > 0) && (y < NY_ht) && (x < NX_ht)))
               { frgb_t clr = fvw_texture_get_color(tx, x, y, (vom+voo)/2, vmin_cur, vmax_cur);
                 fvw_paint_vert_square(xlo, ylo, xhi, ylo, vom, voo, zscale, clr.c[0], clr.c[1], clr.c[2]);
               }
	     /* Step to next pixel in same row: */
             vmo = voo; pvo++;
             pvm++;
	   }  
	glEnd();
        /* Swap row buffers, {vm<-->vo}: */
        { float *t = vo; vo = vm; vm = t; }
      }
    if (fvw_debug_paint) { fprintf(stderr, "- %s\n", __FUNCTION__); }
  }
