/* See langev_util.h. */
/* Last edited on 2024-12-21 14:00:06 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

#include <jsrandom.h>
#include <affirm.h>
#include <bool.h>
#include <uint16_image.h>
#include <vec.h>

#include <langev_util.h>

#define TRUNCATED_EXP_BIAS (0.00004539992976248485)
  /* The value {exp(TRUNCATED_EXP_MIN_Z)}. */

double truncated_exp(double z)
  {
    if (z < TRUNCATED_EXP_MIN_Z) 
      { return 0.0; }
    else
      { return exp(z) - TRUNCATED_EXP_BIAS; }
  }

#define TRUNCATED_BELL_BIAS (0.00001525878906250000)
  /* The value {exp(-ln(2)*TRUNCATED_BELL_MAX_Z^2)}. */

double truncated_bell(double z)
  {
    if (fabs(z) > TRUNCATED_BELL_MAX_Z) 
      { return 0.0; }
    else
      { return exp(-M_LN2*z*z) - TRUNCATED_BELL_BIAS; }
  }

#define TRUNCATED_SIGMOID_SCALE (0.9999999845827420852373)
  /* The value {erf(TRUNCATED_SIGMOID_MAX_Z)}. */

double truncated_sigmoid(double z)
  {
    if (fabs(z) > TRUNCATED_SIGMOID_MAX_Z) 
      { return 0.0; }
    else
      { return (1 + erf(z)/TRUNCATED_SIGMOID_SCALE)/2; }
  }

int pick_elem(double w[], int n, int idef, double wtot)
  { 
    /* Paranoia: */
    if (n <= 0) { return idef; }
    int k;
    double wsum; 
    if (wtot > 0)
      { wsum = wtot; }
    else
      { /* Compute the sum {wsum} of the weights: */
        wsum = 0;
        for (k = 0; k < n; k++) { wsum += w[k]; }
        if (wsum <= 0) { return idef; }
      }
    /* Throw the die: */
    double wsel = wsum*drandom();
    wsum = 0;
    for (k = 0; k < n; k++) 
      { double wk = w[k];
        demand(wk >= 0, "negative weight"); 
        wsum += wk; 
        if (wsel <= wsum) { return k; }
      } 
    /* Roundoff errors or a bad {wtot} may bring us here: */
    return idef;
  }

void pick_colors(int M, int32_t x[], int n, float rgb[])
  { /* Create a matrix {mat[3*r+k]}, for {r} in {0..M-1} and {k} in {0..2}, that
      projects the {M}-cube corners {{-1,+1}^ M} into the unit ball of {R^ 3}. */
    double mat[3*M];
    int r, k, i;
    for (r = 0; r < M; r++)
      { mat[3*r+0] = 1.0;
        mat[3*r+1] = cos(2*M_PI*r/M);
        mat[3*r+2] = sin(2*M_PI*r/M);
      } 
    /* Normalize each column so that its length is {1/sqrt(M)}: */
    for (k = 0; k < 3; k++)
      { double sii = 0; 
        for (r = 0; r < M; r++) { double mi = mat[3*r+k]; sii += mi*mi; }
        sii = sqrt(M*sii);
        for (r = 0; r < M; r++) { mat[3*r+k] /= sii; }
      }
    /* Map trait vectors to colors: */
    for (i = 0; i < n; i++)
      { for (k = 0; k < 3; k++)
          { int32_t xi = x[i];
            double cik = 0;
            for (r = 0; r < M; r++)
              { double xir = 2.0*(xi & 1) - 1.0;
                cik += xir*mat[3*r+k];
                xi /= 2;
             }
            rgb[3*i+k] = cik;
          }
      }
  }

uint16_image_t *colorize_frame(frame_t *frm, sibs_color_func_t *color)
  { int nx = frm->cols;
    int ny = frm->rows;
    /* Allocate a color image: */
    int chns = 3;
    uint16_image_t *img = uint16_image_new(nx, ny, chns);
    img->maxval = uint16_image_MAX_SAMPLE;
    int maxval = (int)img->maxval;
    /* Fill it with colors: */
    int y, x;
    for (y = 0; y < ny; y++)
      { uint16_t *smp = img->smp[ny - 1 - y];
        for (x = 0; x < nx; x++)
          { float rgb[3]; /* RGB color, [0_1] scale. */
            /* Get site occupation in {frm}: */
            site_id_t p = site_id(x,y, nx,ny);
            site_t *c = get_site_address(frm, p);
            color(p, &(c->oc), rgb);
            int ich;
            for (ich = 0; ich < chns; ich++)
              { double rgbi = rgb[ich];
                assert(rgbi >= 0.0);
                assert(rgbi <= 1.0);
                int smpi = (rgbi == 1.0 ? maxval : (int)floor(rgbi*(maxval + 1)));
                assert((smpi >= 0) && (smpi <= maxval));
                (*smp) = (uint16_t)smpi;
                smp++;
              }
          }

      }
    return img;
  }

void print_hex_alpha(FILE *wr, uint32_t w, int nbits)
  { if (nbits == 0) { return; }
    int ndigs = (nbits + 3)/4;
    uint32_t d = 1u << (4*(ndigs-1));
    while (d != 0)
      { uint32_t q = w/d;
        fputc('A' + q, wr);
        w %= d;
        d >>= 4;
      }
  }

void print_binary(FILE *wr, uint32_t w, int nbits)
  { if (nbits == 0) { return; }
    gene_t mask = (1u << (nbits-1));
    while(mask != 0)
      { /* Print bit: */
        fputc(((w & mask) != 0 ? '1' : '0'), wr);
        /* Go to next bit: */
        mask >>= 1;
      }
  }
