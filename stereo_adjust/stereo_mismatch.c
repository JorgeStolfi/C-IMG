// See stereo-mismatch.h.
// Last edited on 2023-11-26 06:49:20 by stolfi

#include <stereo_mismatch.h>

#include <bool.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

// INTERNAL PROTOTYPES

void stm_add_12_mismatch_term
  ( bool_t grad, bool_t verbose, 
    r2_t *p1_p, r2_t *p2_p, stereo_params_t *L_p, double zt1, double zt2,
    double *S_p, stereo_params_t *DS_DL_p, double DS_Dzt1_p[], double DS_Dzt2_p[]
  );
/* 
  Adds to the total mismatch {*S_p} the mismatch between the predicted
  position of a point in image 2, and its observed position {*p2_p =
  (xp2,yp2)} and assumed height {zt2}. The prediction uses the
  observed PCS point coordinates {*p1_p = (xp1,yp1)} in image 1, its
  presumed ICS height {zt1}, and the current 1->2 transformation
  parameters {*L_p}. If {grad = TRUE}, the procedure also accumulates
  the derivatives of that term to the derivatives of {*S_p}. */

void stm_add_21_mismatch_term
  ( bool_t grad, bool_t verbose, 
    r2_t *p1_p, r2_t *p2_p, stereo_params_t *L_p, double zt1, double zt2,
    double *S_p, stereo_params_t *DS_DL_p, double DS_Dzt1_p[], double DS_Dzt2_p[]
  );
/*
  Adds to the total mismatch {*S_p} the mismatch between the predicted
  position of a point in image 1, and its actual position{*p1_p =
  (xp1,yp1)} and assumed height {zt1}. The prediction uses The
  prediction uses the observed PCS point coordinates {*p2_p =
  (xp2,yp2)} in image 2, its presumed ICS height {zt2}, and the
  current 1->2 transformation parameters {*L_p}, in the inverse
  direction. If {grad = TRUE}, the procedure also accumulates the
  derivatives of that term to the derivatives of {*S_p}. */
  
void stm_scale_mismatch(bool_t grad, int n, double *S_p, stereo_params_t *DS_DL_p, double DS_Dzt1[], double DS_Dzt2[]);
/*
  Scales the mismatch {*S_p} by the factor {1/n}.
  If {grad = TRUE}, also scales the gradient. */
  
void stm_make_tangential(r3x3_t *DR_p, r3x3_t *R_p);
/*
  Projects each row of {*DR_p} so that it is orthogonal to the 
  corresponding row of {*R_p}. */

void stm_modify_param_abs(double *x_p, double e, double max_ch_abs, double *ch_abs_p);
/*
  Increments variable {*x_p} by {e}, clipping the latter to {±max_ch_abs}.
  Also adds the square of the clipped increment to {*ch_abs_p}. */

void stm_modify_param_rel(double *x_p, double e, double max_ch_rel, double *ch_rel_p);
/*
  Multiplies {*x_p} by a factor {exp(e/(*x_p))}, clipping that factor to
  between {exp(max_ch_rel)} and {exp(-max_ch_rel)}. Note that this is
  approximately equivalent to adding {e} to {*x_p}, for small {e}. The
  procedure also adds the square of the log of the clipped scale
  factor to {*ch_rel_p}. */

void stm_test_grad_comp
  ( int n, r2_t *p1, r2_t *p2, stereo_params_t *L_p, double zt1[], double zt2[],
    double S_a, double h, char *parm_name, double *parm_p, double deriv_grd
  );
/*
  Compares a component {deriv_grad} of the gradient with the
  corresponding numerical derivative, computed by varying the
  parameter {*parm_p} (which must be some component of {*L_p},
  {zt1[k]}, or {zt2[k]}) by {h} and recomputing the mismatch. Assumes
  that the current mismatch is {S_a}. The results are printed using
  {*parm_name} as the component's name. The vale of {*parm_p} is
  restored before exiting. */

// DEFINITIONS

void stm_invert_XYZ_perspective(r3_t *p_p, stereo_params_t *L_p, r3_t *t_p)
{
  double d = L_p->d;
  double f = d/(d + p_p->c[2]);
  t_p->c[0] = p_p->c[0] * f;
  t_p->c[1] = p_p->c[1] * f;
  t_p->c[2] = p_p->c[2] * f;
}

void stm_invert_XY_perspective(r2_t *p_p, double zt, stereo_params_t *L_p, r2_t *t_p)
{
  double d = L_p->d;
  double f = (d - zt)/d;
  t_p->c[0] = p_p->c[0] * f;
  t_p->c[1] = p_p->c[1] * f;
}

void stm_guess_translation(int n, r2_t *p1, r2_t *p2, stereo_params_t *L_p, double zt1[], double zt2[])
{  
  r3_t c1, c2, r1, s1;
  r2_t t1, t2;
  int i,j,k;
  for (j = 0; j < 3; j++) { c1.c[j] = 0; c2.c[j] = 0; }
  
  for (k = 0; k < n; k++)
    { r2_t p1k = (r2_t){{p1[k].c[0], p1[k].c[1]}};
      r2_t p2k = (r2_t){{p2[k].c[0], p2[k].c[1]}};
      stm_invert_XY_perspective(&p1k, zt1[k], L_p, &t1);
      stm_invert_XY_perspective(&p2k, zt2[k], L_p, &t2);
      c1.c[0] += t1.c[0]; c2.c[0] += t2.c[0];
      c1.c[1] += t1.c[1]; c2.c[1] += t2.c[1];
      c1.c[2] += zt1[k];  c2.c[2] += zt2[k];
    }
  
  for (j = 0; j < 3; j++) { c1.c[j] /= n; c2.c[j] /= n; }
  
  // fprintf(stderr, "  R =\n");
  // print_r3x3_t(stderr, "    ", &(L_p->R), "\n");
  
  for (j = 0; j < 3; j++) 
    { r1.c[j] = 0;
      for (i = 0; i < 3; i++)
        { r1.c[j] += c1.c[i] * L_p->R.c[i][j]; }
      s1.c[j] = r1.c[j] * L_p->m;
    }

  fprintf(stderr, "  c1 = ( %6.1f %6.1f %6.1f )\n", c1.c[0], c1.c[1], c1.c[2]);
  fprintf(stderr, "  r1 = ( %6.1f %6.1f %6.1f )\n", r1.c[0], r1.c[1], r1.c[2]);
  fprintf(stderr, "  s1 = ( %6.1f %6.1f %6.1f )\n", s1.c[0], s1.c[1], s1.c[2]);
  fprintf(stderr, "  c2 = ( %6.1f %6.1f %6.1f )\n", c2.c[0], c2.c[1], c2.c[2]);
  fprintf(stderr, "\n");

  for (j = 0; j < 3; j++) { L_p->v.c[j] = c2.c[j] - s1.c[j]; }
}

void stm_compute_mismatch(bool_t grad, bool_t verbose, 
    int n, r2_t *p1, r2_t *p2, stereo_params_t *L_p, double zt1[], double zt2[],
    double *S_p, stereo_params_t *DS_DL_p, double DS_Dzt1[], double DS_Dzt2[])
{
  int i, j, k;

  if (verbose)
    { fprintf(stderr, "\n");
      fprintf(stderr, "--- enter stm_compute_mismatch -------------\n");
      stm_print_corresp(stderr, n, L_p, zt1, zt2);
      fprintf(stderr, "\n");
    }

  (*S_p) = 0;
  if (grad) 
    { // Clear gradient accumulators
      for ( k = 0; k < n; k++) { DS_Dzt1[k] = 0; DS_Dzt2[k] = 0; }
      DS_DL_p->d = 0;
      for (i = 0; i < 3; i++)
        { DS_DL_p->v.c[i] = 0;
          for (j = 0; j < 3; j++)
            { DS_DL_p->R.c[i][j] = 0; }
        }
      DS_DL_p->m = 0;
    }
  for (k = 0; k < n; k++)
    { 
      stm_add_12_mismatch_term(grad, FALSE,
        &(p1[k]), &(p2[k]), L_p, zt1[k], zt2[k], 
        S_p, DS_DL_p, &(DS_Dzt1[k]), &(DS_Dzt2[k]));
      stm_add_21_mismatch_term(grad, FALSE,
        &(p1[k]), &(p2[k]), L_p, zt1[k], zt2[k], 
        S_p, DS_DL_p, &(DS_Dzt1[k]), &(DS_Dzt2[k]));
    }
  stm_scale_mismatch(grad, n, S_p, DS_DL_p, DS_Dzt1, DS_Dzt2);

  if (verbose)
    { fprintf(stderr, "\n");
      fprintf(stderr, "mismatch = %12.5f\n", *S_p);
      if (grad) 
        { fprintf(stderr, "gradient = \n");
          stm_print_corresp(stderr, n, DS_DL_p, DS_Dzt1, DS_Dzt2);
          for (k = 0; k < n; k++)
            { fprintf(stderr, "  %12.5f %12.5f\n", DS_Dzt1[k], DS_Dzt2[k]); }
        }
      fprintf(stderr, "--- exit stm_compute_mismatch --------------\n");
    }
}

void stm_add_12_mismatch_term(bool_t grad, bool_t verbose, 
    r2_t *p1_p, r2_t *p2_p, stereo_params_t *L_p, double zt1, double zt2,
    double *S_p, stereo_params_t *DS_DL_p, double DS_Dzt1_p[], double DS_Dzt2_p[])
{
  double f1, xt1, yt1, xr, yr, zr, xs, ys, zs, xw2, yw2, zw2, f2, xh2, yh2, xe2, ye2, ze2;
  double DS_Df1, DS_Dxt1, DS_Dyt1, DS_Dxr, DS_Dyr, DS_Dzr, DS_Dxs, DS_Dys, DS_Dzs, 
         DS_Dxw2, DS_Dyw2, DS_Dzw2, DS_Df2, DS_Dxh2, DS_Dyh2, DS_Dxe2, DS_Dye2, DS_Dze2;

  double d = L_p->d;
  double m = L_p->m;
  double S = 0;
  double DS_Dd = 0, DS_Dm = 0;
  double DS_Dzt1 = 0, DS_Dzt2 = 0;

  f1 = (d - zt1)/d; DS_Df1 = 0;

  xt1 = p1_p->c[0] * f1; DS_Dxt1 = 0;
  yt1 = p1_p->c[1] * f1; DS_Dyt1 = 0;

  if (verbose) { fprintf(stderr, "  t1 = ( %6.1f %6.1f %6.1f )\n", xt1, yt1, zt1); }

  xr = L_p->R.c[0][0]*xt1 + L_p->R.c[1][0]*yt1 + L_p->R.c[2][0]*zt1; DS_Dxr = 0;
  yr = L_p->R.c[0][1]*xt1 + L_p->R.c[1][1]*yt1 + L_p->R.c[2][1]*zt1; DS_Dyr = 0;
  zr = L_p->R.c[0][2]*xt1 + L_p->R.c[1][2]*yt1 + L_p->R.c[2][2]*zt1; DS_Dzr = 0;

  xs = m*xr; DS_Dxs = 0;
  ys = m*yr; DS_Dys = 0;
  zs = m*zr; DS_Dzs = 0;

  xw2 = L_p->v.c[0] + xr; DS_Dxw2 = 0;
  yw2 = L_p->v.c[1] + yr; DS_Dyw2 = 0;
  zw2 = L_p->v.c[2] + zr; DS_Dzw2 = 0;
  
  if (verbose) { fprintf(stderr, "  w2 = ( %6.1f %6.1f %6.1f )\n", xw2, yw2, zw2); }

  f2 = d/(d - zw2);  DS_Df2 = 0;

  xh2 = xw2 * f2;  DS_Dxh2 = 0;
  yh2 = yw2 * f2;  DS_Dyh2 = 0;

  if (verbose) { fprintf(stderr, "  h2 = ( %6.1f %6.1f )\n", xh2, yh2); }
  if (verbose) { fprintf(stderr, "  p2 = ( %6.1f %6.1f )\n", p2_p->c[0], p2_p->c[1]); }

  xe2 = xh2 - p2_p->c[0];  DS_Dxe2 = 0;
  ye2 = yh2 - p2_p->c[1];  DS_Dye2 = 0;
  ze2 = zw2 - zt2;         DS_Dze2 = 0;

  if (verbose) { fprintf(stderr, "  e2 = ( %6.1f %6.1f %6.1f )\n\n", xe2, ye2, ze2); }

  S += xe2*xe2 + ye2*ye2 + ze2*ze2;
  
  if (grad)
    { DS_Dxe2 += 2*xe2; DS_Dye2 += 2*ye2; DS_Dze2 += 2*ze2;

      DS_Dzw2 += DS_Dze2; DS_Dzt2 += -DS_Dze2;
      DS_Dyh2 += DS_Dye2;
      DS_Dxh2 += DS_Dxe2;

      DS_Dyw2 += DS_Dyh2*f2; DS_Df2 += DS_Dyh2*yw2;
      DS_Dxw2 += DS_Dxh2*f2; DS_Df2 += DS_Dxh2*xw2;

      DS_Dd += DS_Df2*f2*(1 - f2)/d; 
      DS_Dzw2 += DS_Df2*f2*f2/d;

      DS_DL_p->v.c[2] += DS_Dzw2; DS_Dzr += DS_Dzw2;
      DS_DL_p->v.c[1] += DS_Dyw2; DS_Dyr += DS_Dyw2;
      DS_DL_p->v.c[0] += DS_Dxw2; DS_Dxr += DS_Dxw2;

      DS_Dzr += DS_Dzs*m;  DS_Dm += DS_Dzs*zr;
      DS_Dyr += DS_Dys*m;  DS_Dm += DS_Dys*yr;
      DS_Dxr += DS_Dxs*m;  DS_Dm += DS_Dxs*xr;

      DS_DL_p->R.c[0][2] += DS_Dzr*xt1;  DS_Dxt1 += DS_Dzr*L_p->R.c[0][2];
      DS_DL_p->R.c[1][2] += DS_Dzr*yt1;  DS_Dyt1 += DS_Dzr*L_p->R.c[1][2];
      DS_DL_p->R.c[2][2] += DS_Dzr*zt1;  DS_Dzt1 += DS_Dzr*L_p->R.c[2][2];

      DS_DL_p->R.c[0][1] += DS_Dyr*xt1;  DS_Dxt1 += DS_Dyr*L_p->R.c[0][1];
      DS_DL_p->R.c[1][1] += DS_Dyr*yt1;  DS_Dyt1 += DS_Dyr*L_p->R.c[1][1];
      DS_DL_p->R.c[2][1] += DS_Dyr*zt1;  DS_Dzt1 += DS_Dyr*L_p->R.c[2][1];

      DS_DL_p->R.c[0][0] += DS_Dxr*xt1;  DS_Dxt1 += DS_Dxr*L_p->R.c[0][0];
      DS_DL_p->R.c[1][0] += DS_Dxr*yt1;  DS_Dyt1 += DS_Dxr*L_p->R.c[1][0];
      DS_DL_p->R.c[2][0] += DS_Dxr*zt1;  DS_Dzt1 += DS_Dxr*L_p->R.c[2][0];

      DS_Df1 += DS_Dyt1 * p1_p->c[1];
      DS_Df1 += DS_Dxt1 * p1_p->c[0];

      DS_Dd += DS_Df1*(1 - f1)/d;  DS_Dzt1 += -DS_Df1/d;
    
      DS_DL_p->m += DS_Dm;
      DS_DL_p->d += DS_Dd;
      (*DS_Dzt1_p) += DS_Dzt1;
      (*DS_Dzt2_p) += DS_Dzt2;
    }

  (*S_p) += S;
}

void stm_add_21_mismatch_term(bool_t grad, bool_t verbose,
    r2_t *p1_p, r2_t *p2_p, stereo_params_t *L_p, double zt1, double zt2,
    double *S_p, stereo_params_t *DS_DL_p, double DS_Dzt1_p[], double DS_Dzt2_p[])
{
  double f2, xt2, yt2, xr, yr, zr, xs, ys, zs, xw1, yw1, zw1, f1, xh1, yh1, xe1, ye1, ze1;
  double DS_Df2, DS_Dxt2, DS_Dyt2, DS_Dxr, DS_Dyr, DS_Dzr, DS_Dxs, DS_Dys, DS_Dzs,
         DS_Dxw1, DS_Dyw1, DS_Dzw1, DS_Df1, DS_Dxh1, DS_Dyh1, DS_Dxe1, DS_Dye1, DS_Dze1;

  double d = L_p->d;
  double m = L_p->m;
  double S = 0;
  double DS_Dd = 0, DS_Dm = 0;
  double DS_Dzt1 = 0, DS_Dzt2 = 0;

  f2 = (d - zt2)/d;
  DS_Df2 = 0;

  xt2 = p2_p->c[0] * f2; DS_Dxt2 = 0;
  yt2 = p2_p->c[1] * f2; DS_Dyt2 = 0;

  if (verbose) { fprintf(stderr, "  t2 = ( %6.1f %6.1f %6.1f )\n", xt2, yt2, zt2); }

  xs = xt2 - L_p->v.c[0]; DS_Dxs = 0;
  ys = yt2 - L_p->v.c[1]; DS_Dys = 0;
  zs = zt2 - L_p->v.c[2]; DS_Dzs = 0;

  xr = xs/m; DS_Dxr = 0;
  yr = ys/m; DS_Dyr = 0;
  zr = zs/m; DS_Dzr = 0;

  xw1 = L_p->R.c[0][0]*xr + L_p->R.c[0][1]*yr + L_p->R.c[0][2]*zr; DS_Dxw1 = 0;
  yw1 = L_p->R.c[1][0]*xr + L_p->R.c[1][1]*yr + L_p->R.c[1][2]*zr; DS_Dyw1 = 0;
  zw1 = L_p->R.c[2][0]*xr + L_p->R.c[2][1]*yr + L_p->R.c[2][2]*zr; DS_Dzw1 = 0;

  if (verbose) { fprintf(stderr, "  w1 = ( %6.1f %6.1f %6.1f )\n", xw1, yw1, zw1); }

  f1 = d/(d - zw1);  DS_Df1 = 0;

  xh1 = xw1 * f1;  DS_Dxh1 = 0;
  yh1 = yw1 * f1;  DS_Dyh1 = 0;

  if (verbose) { fprintf(stderr, "  h1 = ( %6.1f %6.1f )\n", xh1, yh1); }
  if (verbose) { fprintf(stderr, "  p1 = ( %6.1f %6.1f )\n", p1_p->c[0], p1_p->c[1]); }

  xe1 = xh1 - p1_p->c[0];  DS_Dxe1 = 0;
  ye1 = yh1 - p1_p->c[1];  DS_Dye1 = 0;
  ze1 = zw1 - zt1;         DS_Dze1 = 0;

  if (verbose) { fprintf(stderr, "  e1 = ( %6.1f %6.1f %6.1f )\n\n", xe1, ye1, ze1); }

  S += xe1*xe1 + ye1*ye1 + ze1*ze1;
  
  if (grad)
    { 
      DS_Dxe1 += 2*xe1; DS_Dye1 += 2*ye1; DS_Dze1 += 2*ze1;

      DS_Dzw1 += DS_Dze1; DS_Dzt1 += -DS_Dze1;
      DS_Dyh1 += DS_Dye1;
      DS_Dxh1 += DS_Dxe1;

      DS_Dyw1 += DS_Dyh1*f1; DS_Df1 += DS_Dyh1*yw1;
      DS_Dxw1 += DS_Dxh1*f1; DS_Df1 += DS_Dxh1*xw1;

      DS_Dd += DS_Df1*f1*(1 - f1)/d; 
      DS_Dzw1 += DS_Df1*f1*f1/d;

      DS_DL_p->R.c[2][0] += DS_Dzw1*xr;  DS_Dxr += DS_Dzw1*L_p->R.c[2][0];
      DS_DL_p->R.c[2][1] += DS_Dzw1*yr;  DS_Dyr += DS_Dzw1*L_p->R.c[2][1];
      DS_DL_p->R.c[2][2] += DS_Dzw1*zr;  DS_Dzr += DS_Dzw1*L_p->R.c[2][2];

      DS_DL_p->R.c[1][0] += DS_Dyw1*xr;  DS_Dxr += DS_Dyw1*L_p->R.c[1][0];
      DS_DL_p->R.c[1][1] += DS_Dyw1*yr;  DS_Dyr += DS_Dyw1*L_p->R.c[1][1];
      DS_DL_p->R.c[1][2] += DS_Dyw1*zr;  DS_Dzr += DS_Dyw1*L_p->R.c[1][2];

      DS_DL_p->R.c[0][0] += DS_Dxw1*xr;  DS_Dxr += DS_Dxw1*L_p->R.c[0][0];
      DS_DL_p->R.c[0][1] += DS_Dxw1*yr;  DS_Dyr += DS_Dxw1*L_p->R.c[0][1];
      DS_DL_p->R.c[0][2] += DS_Dxw1*zr;  DS_Dzr += DS_Dxw1*L_p->R.c[0][2];

      DS_Dzs += DS_Dzr/m;  DS_Dm += -DS_Dzr*zr/m;
      DS_Dys += DS_Dyr/m;  DS_Dm += -DS_Dyr*yr/m;
      DS_Dxs += DS_Dxr/m;  DS_Dm += -DS_Dxr*xr/m;

      DS_DL_p->v.c[2] += -DS_Dzs;  DS_Dzt2 += DS_Dzs;
      DS_DL_p->v.c[1] += -DS_Dys;  DS_Dyt2 += DS_Dys;
      DS_DL_p->v.c[0] += -DS_Dxs;  DS_Dxt2 += DS_Dxs;

      DS_Df2 += DS_Dyt2 * p2_p->c[1];
      DS_Df2 += DS_Dxt2 * p2_p->c[0];

      DS_Dd += DS_Df2*(1 - f2)/d;  DS_Dzt2 += -DS_Df2/d;

      DS_DL_p->m += DS_Dm;
      DS_DL_p->d += DS_Dd;
      (*DS_Dzt1_p) += DS_Dzt1;
      (*DS_Dzt2_p) += DS_Dzt2;
   }

  (*S_p) += S;
}

void stm_scale_mismatch(bool_t grad, int n, double *S_p, stereo_params_t *DS_DL_p, double DS_Dzt1[], double DS_Dzt2[])
{ 
  int i,j,k;
  double s = 1.0/((double)n);
  (*S_p) *= s;
  
  if (grad)
    { 
      DS_DL_p->d *= s;
      DS_DL_p->m *= s;
      for (i = 0; i < 3; i++)
        { DS_DL_p->v.c[i] *= s;
          for (j = 0; j < 3; j++)
            { DS_DL_p->R.c[i][j] *= s; }
        }

      for (k = 0; k < n; k++)
        { DS_Dzt1[k] *= s; DS_Dzt2[k] *= s; }
    }
}

double stm_gradient_dot(int n, 
    stereo_params_t *UL_p, double Uzt1[], double Uzt2[], 
    stereo_params_t *VL_p, double Vzt1[], double Vzt2[])
{
  double x, y;
  double s = 0;
  int i, j, k;
  
  x = UL_p->d; y = VL_p->d; s += x*y;
  for (i = 0; i < 3; i++) 
    { for (j = 0; j < 3; j++) 
        { x = UL_p->R.c[i][j]; y = VL_p->R.c[i][j]; s += x*y; }
    }
  x = UL_p->m; y = VL_p->m; s += x*y;
  for (j = 0; j < 3; j++) 
    { x = UL_p->v.c[j]; y = VL_p->v.c[j]; s += x*y; }
  for (k = 0; k < n; k++)
    { x = Uzt1[k]; y = Vzt1[k]; s += x*y;
      x = Uzt2[k]; y = Vzt2[k]; s += x*y;
    }
  return s;
}

double stm_gradient_norm_sqr(int n, stereo_params_t *UL_p, double Uzt1[], double Uzt2[])
{
  double x;
  double s = 0;
  int i, j, k;
  
  x = UL_p->d; s += x*x;
  for (i = 0; i < 3; i++) 
    { for (j = 0; j < 3; j++) 
        { x = UL_p->R.c[i][j]; s += x*x; }
    }
  x = UL_p->m; s += x*x;
  for (j = 0; j < 3; j++) 
    { x = UL_p->v.c[j]; s += x*x; }
  for (k = 0; k < n; k++)
    { x = Uzt1[k]; s += x*x;
      x = Uzt2[k]; s += x*x;
    }
  return s;
}

void stm_modify_params(double h, 
  double *ch_rel_p, double *ch_abs_p, 
  double max_ch_rel, double max_ch_abs,
  int n, stereo_params_t *L_p, double zt1[], double zt2[],
  stereo_params_t *DS_DL_p, double DS_Dzt1[], double DS_Dzt2[])
{
  int i, j, k;
  (*ch_abs_p) = 0; 
  (*ch_rel_p) = 0;
  for (k = 0; k < n; k++)
    { stm_modify_param_abs(&(zt1[k]), -h*DS_Dzt1[k], max_ch_abs, ch_abs_p);
      stm_modify_param_abs(&(zt2[k]), -h*DS_Dzt2[k], max_ch_abs, ch_abs_p);
    }
  stm_modify_param_rel(&(L_p->d), -h*DS_DL_p->d, max_ch_rel, ch_rel_p);
  // Should turn the gradient into a rotation by exponentiating it.
  // Hack: make the gradient tangential to {R}, add and orthonormalize.
  // However, note that the change must be added to {*ch_rel_p} not {*ch_abs_p}
  stm_make_tangential(&(DS_DL_p->R), &(L_p->R));
  for (i = 0; i < 3; i++)
    { for (j = 0; j < 3; j++)
        { stm_modify_param_abs(&(L_p->R.c[i][j]), -h*DS_DL_p->R.c[i][j], max_ch_rel, ch_rel_p); }
    }
  orthonormalize(&(L_p->R));
  stm_modify_param_rel(&(L_p->m), -h*DS_DL_p->m, max_ch_rel, ch_rel_p);
  for (i = 0; i < 3; i++)
    { stm_modify_param_abs(&(L_p->v.c[i]), -h*DS_DL_p->v.c[i], max_ch_abs, ch_abs_p); }
  (*ch_abs_p) = sqrt((*ch_abs_p));
  (*ch_rel_p) = sqrt((*ch_rel_p));
}

void stm_modify_param_abs(double *x_p, double e, double max_ch_abs, double *ch_abs_p)
{  
  if (e < -max_ch_abs) { e = -max_ch_abs; }
  if (e > max_ch_abs) { e = max_ch_abs; }
  (*x_p) += e; (*ch_abs_p) += e*e;
}

void stm_modify_param_rel(double *x_p, double e, double max_ch_rel, double *ch_rel_p)
{  
  e /= (*x_p);
  if (e < -max_ch_rel) { e = -max_ch_rel; }
  if (e > max_ch_rel) { e = max_ch_rel; }
  (*x_p) *= exp(e); (*ch_rel_p) += e*e;
}

void stm_make_tangential(r3x3_t *DR_p, r3x3_t *R_p)
{  
  int i, j;
  double s;
  for (i = 0; i < 3; i++)
    { // Orthogonalize row i against preceding rows
      s = 0; 
      for (j = 0; j < 3; j++) { s += DR_p->c[i][j] * R_p->c[i][j]; }
      for (j = 0; j < 3; j++) { DR_p->c[i][j] -= s*R_p->c[i][j]; }
  }
}

void stm_print_params(FILE *f, stereo_params_t *L_p)
{
  fprintf(f, "  d = %6.1f\n", L_p->d);
  fprintf(f, "  R = \n");
  print_r3x3_t(f, "    ", &(L_p->R), "\n");
  fprintf(f, "  m = %8.5f\n", L_p->m);
  fprintf(f, "  v = ");
  print_r3_t(stderr, "", &(L_p->v), "\n");
}

void stm_print_corresp(FILE *f, int n, stereo_params_t *L_p, double zt1[], double zt2[])
{
  int k;
  stm_print_params(f, L_p);
  for (k = 0; k < n; k++)
    { fprintf(stderr, "  %12.5f %12.5f\n", zt1[k], zt2[k]); }
}

void stm_test_gradient(int n, r2_t *p1, r2_t *p2, stereo_params_t *L_p, double zt1[], double zt2[])
{
  double h = 1.0/1024.0/1024.0;
  char lab[20];

  double  S_a;
  stereo_params_t DS_DL;
  double *DS_Dzt1 = talloc(n, double);
  double *DS_Dzt2 = talloc(n, double);
  
  int i,j,k;
  
  fprintf(stderr, "\n=== BEGIN GRADIENT TEST =====================\n");

  stm_compute_mismatch(TRUE, FALSE, n, p1, p2, L_p, zt1, zt2, &S_a, &DS_DL, DS_Dzt1, DS_Dzt2);

  stm_test_grad_comp(n, p1, p2, L_p, zt1, zt2, S_a, h, "d", &(L_p->d), DS_DL.d);
  for (i = 0; i < 3; i++)
    { for (j = 0; j < 3; j++) 
        { sprintf(lab, "R[%d,%d]", i,j);
          stm_test_grad_comp(n, p1, p2, L_p, zt1, zt2, S_a, h, lab, &(L_p->R.c[i][j]), DS_DL.R.c[i][j]);
        }
    }
  stm_test_grad_comp(n, p1, p2, L_p, zt1, zt2, S_a, h, "m", &(L_p->m), DS_DL.m);
  for (i = 0; i < 3; i++)
    { sprintf(lab, "v[%d]", i);
      stm_test_grad_comp(n, p1, p2, L_p, zt1, zt2, S_a, h, lab, &(L_p->v.c[i]), DS_DL.v.c[i]);
    }
  for (k = 0; k < n; k++)
    { sprintf(lab, "zt1[%d]", k);
      stm_test_grad_comp(n, p1, p2, L_p, zt1, zt2, S_a, h, lab, &(zt1[k]), DS_Dzt1[k]);
      sprintf(lab, "zt2[%d]", k);
      stm_test_grad_comp(n, p1, p2, L_p, zt1, zt2, S_a, h, lab, &(zt2[k]), DS_Dzt2[k]);
    }
  fprintf(stderr, "\n=== END GRADIENT TEST =======================\n");
}

void stm_test_grad_comp(
   int n, r2_t *p1, r2_t *p2, stereo_params_t *L_p, double zt1[], double zt2[],
   double S_a, double h, char *parm_name, double *parm_p, double deriv_grd)
{
  double S_b, deriv_num;
  double old_parm = (*parm_p);
  (*parm_p) += h;
  stm_compute_mismatch(FALSE, FALSE, n, p1, p2, L_p, zt1, zt2, &S_b, NULL, NULL, NULL);
  deriv_num = (S_b - S_a)/h;
  fprintf(stderr, "%-20s num = %12.5f  grd = %12.5f\n", parm_name, deriv_num, deriv_grd);
  (*parm_p) = old_parm;
}

