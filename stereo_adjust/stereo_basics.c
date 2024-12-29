/*  See {stero_basics.h} */
/*  Last edited on 2023-11-26 06:48:44 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <r3.h>
#include <r3x3.h>
#include <stereo_basics.h>

void orthonormalize(r3x3_t *R_p)
{
  double s,v;
  for (uint32_t i = 0;  i < 3; i++)
    { /*  Orthogonalize row {i} against preceding rows */
      for (uint32_t k = 0;  k < i; k++)
        { s = 0; 
          for (uint32_t j = 0;  j < 3; j++) { s += R_p->c[i][j]*R_p->c[k][j]; }
          for (uint32_t j = 0;  j < 3; j++) { R_p->c[i][j] -= s*R_p->c[k][j]; }
        }
      /*  Normalize row {i}: */
      s = 0; 
      for (uint32_t j = 0;  j < 3; j++) { v = R_p->c[i][j]; s += v*v; }
      demand(s > 0.0, "*** degenerate rotation ***\n");
      s = 1.0/sqrt(s);
      for (uint32_t j = 0;  j < 3; j++) { R_p->c[i][j] *= s; }
  }
}

void print_r3_t(FILE *f, char *pref, r3_t *v_p, char *suff)
{
  fprintf(f, "%s(", pref);
  for (uint32_t j = 0;  j < 3; j++)
    { fprintf(f, " %6.1f", v_p->c[j]); }
  fprintf(f, " )%s", suff);
}

void print_r3x3_t(FILE *f, char *pref, r3x3_t *R_p, char *suff)
{
  for (uint32_t i = 0;  i < 3; i++)
    { fprintf(f, "%s(", pref);
      for (uint32_t j = 0;  j < 3; j++)
        { fprintf(f, " %8.5f", R_p->c[i][j]); }
      fprintf(f, " )%s", suff);
    }
}




