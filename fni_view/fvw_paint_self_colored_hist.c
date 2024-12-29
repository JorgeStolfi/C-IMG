/* See fvw_paint_self_colored.h */
/* Last edited on 2024-12-23 09:06:36 by stolfi */

#include <assert.h>
#include <math.h>
#include <GL/glu.h>

#include <float_image.h>
#include <frgb.h>
#include <frgb_path.h>

#include <fvw_paint.h>
#include <fvw_paint_self_colored.h>
#include <fvw_paint_self_colored_hist.h>


void fvw_paint_self_colored_hist_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double zscale, 
    float vmin,
    float vmax,
    bool_t skirts
  )
  {
    if (fvw_debug_paint) { fprintf(stderr, "+ %s\n", __FUNCTION__); }
    
    /* Get  the height image dimensons: */
    assert(ht != NULL);
    int32_t HNC, HNX, HNY;
    float_image_get_size(ht, &HNC, &HNX, &HNY);

    /* We need buffers for two rows of height map samples: */
    float va[HNX], vb[HNX];
    float *vo = va; /* Current row of height samples. */
    float *vm = vb; /* Next row of height samples. */
    
    /* Scan rows of height array: */
    int32_t x, y;
    for(y = 0; y <= HNY; y++)
      { /* If {y > 0}, then {vm} contains the pixel values of row {y-1}. */
	if (y <= HNY-1)
          { /* Get in {vo} the pixel values of row {y}: */
            float_image_get_sample_row(ht, (int32_t)c, 0, HNX-1, y, vo);
          }
	/* Paint the pixels in row {y} and the walls between row {y} and row {y-1}: */
	glBegin(GL_QUADS);
	float *pvo = vo; /* Pointer to pixel value in col {x}, row {y}. */
	float vmo = 0;   /* Pixel value in col {x-1}, row {y} */
        float *pvm = vm; /* Pointer to pixel value in col {x}, row {y-1}. */
	for(x = 0; x <= HNX; x++)
	   { float xlo = (float)x, xhi = (float)(x + 1.0); /* X range of rod. */
             float ylo = (float)y, yhi = (float)(y + 1.0); /* Y range of rof */
             float voo = (float)(x < HNX ? (*pvo) : 0); /* Pixel value in col {x}, row {y} */
             float vom = (float)(x < HNX ? (*pvm) : 0); /* Pixel value in col {x}, row {y-1} */
             /* Paint top of rod: */
	     if(( x < HNX) && ( y < HNY))
	       {fvw_paint_self_colored_horz_square(xlo, xhi, ylo, yhi, voo, zscale, vmin, vmax);}
             /* Paint west wall of rod: */
             if ((y < HNY) && (skirts || ((x > 0) && (x < HNX))))
               { fvw_paint_self_colored_vert_square(xlo, ylo, xlo, yhi, vmo, voo, zscale, vmin, vmax); }
             /* Paint south wall of rod: */
             if (skirts || ((y > 0) && (y < HNY) && (x < HNX)))
               { fvw_paint_self_colored_vert_square(xlo, ylo, xhi, ylo, vom, voo, zscale, vmin, vmax); }
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
  
void fvw_paint_self_colored_horz_square
  ( float xlo, float xhi, 
    float ylo, float yhi, 
    float v,
    double zscale,
    float vmin,
    float vmax 
  )
  {
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
    fvw_color_from_value(3, v - vzer, vdel, clr);
    GLfloat CR = clr[0], CG = clr[1], CB = clr[2];
    glColor3f(CR,CG,CB);
    GLfloat z = (GLfloat)(zscale*v);
    fvw_paint_quad(xlo, ylo, z, xhi, ylo, z, xhi, yhi, z, xlo, yhi, z);
  }

void fvw_paint_self_colored_vert_square
  ( float xa, float ya, 
    float xb, float yb, 
    float vlo, float vhi,
    double zscale,
    float vmin,
    float vmax 
  )
  {
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
    
    float vm = (vlo + vhi)/2;
    float clr[3];
    fvw_color_from_value(3, vm - vzer, vdel, clr);
    GLfloat CR = clr[0], CG = clr[1], CB = clr[2];
    glColor3f(CR,CG,CB);
    GLfloat zlo = (GLfloat)(zscale*vlo);
    GLfloat zhi = (GLfloat)(zscale*vhi);
    fvw_paint_quad(xa, ya, zlo, xb, yb, zlo, xb, yb, zhi, xa, ya, zhi);
  }
