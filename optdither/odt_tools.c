/* See odt_tools.h */
/* Last edited on 2023-03-17 22:07:44 by stolfi */

#define odt_tools_C_COPYRIGHT \
  "Copyright © 1993 by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>

#include <uint16_image.h>
#include <float_image.h>
#include <uint16_image_check_dither.h>
#include <uint16_image_read_gen.h>

#include <odt_tools.h>

float_image_t *odt_read_dither_matrix(char *fname)
  { 
    /* Read the dither matrix as integer valued image {xp}: */
    image_file_format_t ffmt = image_file_format_PNM;
    uint32_t imaxval[uint16_image_read_gen_MAX_CHANS];
    double gamma, bias;
    bool_t verbose_read = TRUE;
    uint16_image_t *xp = uint16_image_read_gen_named
      ( fname, ffmt, imaxval, &gamma, &bias, verbose_read );
    demand(xp->chns == 1, "dither matrix must be grayscale");
    assert(xp->vaxval == imaxval[0]); 
    if ((gamma != 0) || (bias != 0))
      { /* Dither matrices are asssumed to be in linear lum scale: */
        fprintf(stderr, "nonzero {gamma} or {bias}, ignored"); 
        gamma = 0; bias = 0;
      }
    uint16_t maxval = xp->maxval;
      
    demand(uint16_image_check_dither(xp, FALSE), "not a dither matrix");
      
    bool_t isMask = FALSE; /* Map 0 and {maxval} to midranges. */
    double *lo = NULL; /* Assume vector of zeros. */
    double *hi = NULL; /* Assume vector of ones. */
    bool_t yup = TRUE;
    bool_t verbose_conv = FALSE;
    float_image_t *xf = float_image_from_uint16_image
      ( xp, isMask, lo, hi, yup, verboze_conv);
    
    uint16_image_free(xp);
    (*maxval_P) = maxval;
    return xf;
  }    


/* TO REWRITE */

void odt_round_image(int32_t m, int32_t n, cmp *xc, int32_t *xi)
  { int32_t i,j;
    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        { xi[i*n+j] = (int32_t)floor(xc[i*n+j].re + 0.5); }
  }

void PrT(int32_t m, int32_t n, int32_t x[])
  {
    int32_t i,j;

    for (i = 0; i < m; i++)
      {
        printf("\n");
        for (j = 0; j < n; j++)
          printf(" %3d", x[i*n+j]);
      }
  }

void SaveCT(char *name, int32_t m, int32_t n, cmp *x, char mode, int32_t length)
  {
    char *rname;
    cmp cmpz;
    plr plrz;
    int32_t i;

    rname = calloc(32,sizeof(char));
    strcat(rname,name);
    strcat(rname,".CAdat");
    FILE *f = open_write(rname, TRUE);
    free(rname);

    fprintf(f,"%i %i\n",m,n);

    for(i=0;i<m*n;i++)
      {
        cmpz=*(x+i);
        if (mode=='c')
          fprintf(f,"%f\n%f\n",cmpz.re,cmpz.im);
        else { plrz=mCmp(cmpz);
          fprintf(f,"%f\n%f\n",plrz.r,plrz.an); };
      } /* End for */       

    fclose(f);

    reprint(name,length);

  }

void Dith(int32_t m, int32_t n, float *x, int32_t *y, int32_t mm, int32_t nn, int32_t *matr)
  {
    float dd = 1.0 / mm / nn;

    int32_t ii = 0;
    int32_t jj = 0;
    for (int32_t i = 0; i < m; i++)
      { for (int32_t j = 0; j < n; j++)
          { if (ii == mm)
              { ii = 0; }
            if (jj == nn)
              { jj = 0; }
            if (*(x + i * n + j) > dd * (*(matr + ii * nn + jj) + 0.5))
              { *(y + i * n + j) = 1; }
            else
              { *(y + i * n + j) = 0; }
            jj++;
          };
        ii++;
      }
  }

void rebuild(int32_t m, int32_t n, cmp * x, cmp * y)
  {
    for (int32_t i = 0; i < m / 2; i++)
      for (int32_t j = 0; j < n / 2; j++)
        {
          *(y + i * n + j) = *(x + (i + m / 2) * n + j + n / 2);
          *(y + (i + m / 2) * n + j + n / 2) = *(x + i * n + j);
          *(y + (i + m / 2) * n + j) = *(x + i * n + j + n / 2);
          *(y + i * n + j + n / 2) = *(x + (i + m / 2) * n + j);
        }
  }

void LiveAgua (int32_t m, int32_t n, cmp * x, int32_t *y)
  {
    int32_t *oo = (int32_t *) calloc (m * n, sizeof (int32_t));
    for (int32_t i = 0; i < m; i++)
      for (int32_t j = 0; j < n; j++)
        *(oo + i * n + j) = 0;

    for (int32_t k = m * n; k > 0; k--)
      { float mx = -20000.0;
        for (int32_t i = 0; i < m; i++)
          for (int32_t j = 0; j < n; j++)
            {
              cmp z = *(x + i * n + j);
              if (*(oo + i * n + j) == 0)
                if (z.re >= mx)
                  { mx = z.re;
                    im = i;
                    jm = j;
                  };
            };
        *(y + im * n + jm) = k - 1;
        *(oo + im * n + jm) = 1000;
      };
  }

void MskFlt(int32_t m, int32_t n, cmp * ax, cmp * ay, cmp * az)
  {
    for (int32_t i = 0; i < m; i++)
      for (int32_t j = 0; j < n; j++)
        *(ay + i * n + j) = M (*(ax + i * n + j), *(az + i * n + j));

  }

void Justifier (int32_t m, int32_t n, cmp * x, cmp * y, cmp * res)
  {
    float jst = Power (m, n, x) / Power (m, n, y);
    for (int32_t i = 0; i < m; i++)
      for (int32_t j = 0; j < n; j++)
        {	cmp z = *(y + i * n + j);
          z.re = z.re * jst;
          z.im = z.im * jst;
          *(res + i * n + j) = z;
        }
  }
