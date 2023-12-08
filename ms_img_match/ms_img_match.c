#define PROG_NAME "ms_img_match"
#define PROG_DESC "find camera correspondence matrix and depth between two images"
#define PROG_VERS "1.0"

/* Last edited on 2023-02-25 16:04:45 by stolfi */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  -imageA {IMG[0]} { {X[0][i]} {Y[0][i]} }.. \\\n" \
  "  -imageB {IMG[1]} { {X[1][i]} {Y[1][i]} }.. \\\n" \
  "  [ -sameCamera ] \\\n" \
  "  -outPrefix {OUT_PREFIX} ]" 


#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  This program reads two photo files {IMG[0]} and {IMG[1]} of the same static scene, whose fields of view overlap. The images may have been taken with the same camera or different cameras, from two different positions.  It assumes that each photo is obtained by conical projection.  It estimates the depth of each pixel in each image, and a 3D projective matrix that takes one camera to the other.  It uses a multi-scale image matching approach.\n" \
  "\n" \
  "  The program assumes that the optical axis of each photo is perpendicular to the image plane and passes through its geometric center.  It assumes that the {X} and {Y} coordinates of points on the image are emasured in pixels, and that the scene at each pixel has a {Z} coordinate, also in pixels.  It assumes that the photo was taken with a pinhole camera, with the aperture on the optical axis at some unknown distance {D[k]} ({k} being 0 or 1) above the image plane (also in pixels).  It assumes that each image was produced by appliying an isometry {S[k]} to the scene, then a uniform scaling matrix {U[k]} to map world units to pixels, then a conical perspective distorion {P[k]}, and finally the {Z} coordinate of each pixel was discarded.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -imageA {IMG[0]} { {X[0][i]} {Y[0][i]} }..\n" \
  "  -imageB {IMG[1]} { {X[1][i]} {Y[1][i]} }..\n" \
  "    These arguments specify the two image files to match, and zero or more pairs corresponding points, one in in each image.  The number of points must be the same.  They are treated as hints rather than hard data.\n" \
  "\n" \
  "  -sameCamera}\n" \
  "    This optional argument, if present, specifies that the two images were taking with the same camera, with the same zoom setting, and have not been cropped.\n" \
  "\n" \
  "  -outPrefix {OUT_FILE}\n" \
  "    This mandatory argument specifies the common prefix for all output file names.\n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  whatever(1)\n" \
  "\n" \
  "AUTHOR\n" \
  "  This program was created on 2015-06-14 by J. Stolfi.\n" \
  "MODIFICATION HISTORY\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " dna_seq_filter_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS  
 
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <bool.h>
#include <r2.h>
#include <r3.h>
#include <r3x3.h>
#include <affirm.h>
#include <rn.h>

#include <stereo_basics.h>
#include <stereo_mismatch.h>

typedef struct mim_options_t 
  {  
    struct mim_image_options_t *image[2];  /* Parameters of the two images. */
    bool_t sameCamera;                     /* If true, assumes same camera and zoom. */
    char *outPrefix;
  } mim_options_t;

typedef struct mim_image_options_t 
  {  
    char *fileName;       /* Image file name. */
    r2_t *ref;            /* Reference points. */
  } mim_image_options_t;


/*  PROTOTYPES */

int main(int argc, char **argv);
options_t *mim_parse_options(int argc, char **argv);
mim_image_options_t *mim_parse_image_options(argparser_t *pp);

  
/*  IMPLEMENTATIONS */

int main(int argc, char **argv)
{
  mim_options_t *o = mim_parse_options(argc,argv);
  
  /* Build multiscale pyramids of the two images:*/
  
  
  /* Find starting comparable levels: */
  
  /* Descend: */
  
  stereo_params_t L;
  r2_t *p1, *p2;
  double *zt1;
  double *zt2;
  int n;        /*  Number of point pairs */
  int szx, szy; /*  Image dimensions */
  int i, j, k;

  scanf("%d %d", &szx, &szy);
  
  scanf("%d", &n);

  fprintf(stderr, "sizeof(R) = %d\n", (int)sizeof(L.R));
  
  p1 = (r2_t *)notnull(malloc(n*sizeof(r2_t)), "no mem");
  zt1 = rn_alloc(n);
  p2 = (r2_t *)notnull(malloc(n*sizeof(r2_t)), "no mem");
  zt2 = rn_alloc(n);

  L.d = 300;
  for (i = 0; i < 3; i++)
    { for (j = 0; j < 3; j++)
        { L.R.c[i][j] = (i == j ? 1.0 : 0.0);  }
    }
  L.m = 1;

  for(k = 0; k < n; k++)
    { double x1, y1, x2, y2;
      scanf(" ( %lf %lf ) = ( %lf %lf )", &x1, &y1, &x2, &y2);
      fprintf(stderr, " ( %6.1f %6.1f ) = ( %6.1f %6.1f )\n", x1,y1,x2,y2);
      p1[k] = (r2_t){{x1 - szx/2, szy/2 - y1}}; 
      p2[k] = (r2_t){{x2 - szx/2, szy/2 - y2}};
      zt1[k] = 30*cos(M_PI*k); zt2[k] = 30*cos(M_PI*k);
    }
    
  stm_guess_translation(n, p1, p2, &L, zt1, zt2);
  
  fprintf(stderr, "initial parameters:\n");
  stm_print_params(stderr, &L);
  
  stm_test_gradient(n, p1, p2, &L, zt1, zt2);
 
  adjust_params(n, p1, p2, &L, zt1, zt2);
  
  fprintf(stderr, "final parameters:\n");
  stm_print_params(stderr, &L);
  
  for(k = 0; k < n; k++)
    { 
      double f1, f2, xt1k, yt1k, xt2k, yt2k;
      f1 = (L.d - zt1[k])/L.d;
      xt1k = p1[k].c[0] * f1;
      yt1k = p1[k].c[1] * f1;
 
      f2 = (L.d - zt2[k])/L.d;
      xt2k = p2[k].c[0] * f2;
      yt2k = p2[k].c[1] * f2;

      printf(" ( %6.1f %6.1f %6.1f ) = ( %6.1f %6.1f %6.1f )\n", 
        xt1k, yt1k, zt1[k], xt2k, yt2k, zt2[k]);
    }

  return 0;
}

void adjust_params(int n, r2_t *p1, r2_t *p2, stereo_params_t *L_p, double zt1[], double zt2[])
{ 
  /* Given n pairs of corresponding points {p1[i], p2[i]}, each on one
    image (in pixels), finds the parameters of a perspective
    correspondence that maps points of image 1 to points of image 2. In
    the process, the procedure also determines the Z coordinates (in
    pixels) {zt1[i]} and {zt2[i]} of each point, relative to the
    corresponding image plane, before the perspective transformation.
    Assumes that {*L_p}, {*zt1[i]}, and {*zt2[i]} already contain
    initial guesses for the unknown parameters. */

  /*  Mismatch and gradient at current state: */
  double  S;
  stereo_params_t DS_DL;
  double *DS_Dzt1 = rn_alloc(n);
  double *DS_Dzt2 = rn_alloc(n);
  
  /*  Minimization direction: */
  stereo_params_t UL;
  double *Uzt1 = rn_alloc(n);
  double *Uzt2 = rn_alloc(n);
  
  /*  Tentative second probe state: */
  stereo_params_t L_b;
  double *zt1_b = rn_alloc(n);
  double *zt2_b = rn_alloc(n);
  
  /*  Mismatch and gradient at tentative second state: */
  double S_b;
  stereo_params_t DS_DL_b;
  double *DS_Dzt1_b = rn_alloc(n);
  double *DS_Dzt2_b = rn_alloc(n);
  /*  double *t; */

  /*  Parameter change in last descending step: */
  double ch_abs = 10000.0; double ch_rel = 10000.0;
  
  /*  Maximum parameter change in each component: */
  double max_ch_abs = 300; double max_ch_rel = log(2);

  auto void pack_params(int n, double x[]);
  auto double goal_f(int n, double x[]);
  auto double goal_f(int n, double x[]);
  
  double goal_f(int n, double x[])
    {
      double S;
      stm_compute_mismatch(TRUE, TRUE, n, p1, p2, L_p, zt1, zt2, &S, &DS_DL_b, DS_Dzt1_b, DS_Dzt2_b);
      return S;
    }

  while ((ch_rel > 0.001) || (ch_abs > 0.01))
    { /*  stm_test_gradient(n, p1, p2, L_p, zt1, zt2); */
      stm_compute_mismatch(TRUE, TRUE, n, p1, p2, L_p, zt1, zt2, &S, &DS_DL, DS_Dzt1, DS_Dzt2);
      select_direction(n, p1, p2, L_p, zt1, zt2, &DS_DL, DS_Dzt1, DS_Dzt2, &UL, Uzt1, Uzt2);
      fprintf(stderr, "direction = \n");
      stm_print_corresp(stderr, n, &UL, Uzt1, Uzt2);
      { /*  Try to find the valley bottom in the gradient direction: */
        double b = stm_gradient_dot(n, &DS_DL, DS_Dzt1, DS_Dzt2, &UL, Uzt1, Uzt2);
        double h = 0.1 * 2.0 * S/b;
        int tries = 0;
        int max_tries = 20;
        int k;
        S_b = 2.0 * S + 1.0;
        while (S_b >= S)
          { fprintf(stderr, "modifying: h = %12.8f", h);
            L_b = (*L_p); 
            for (k = 0; k < n; k++) { zt1_b[k] = zt1[k]; zt2_b[k] = zt2[k]; }
            stm_modify_params(h, &ch_rel, &ch_abs, max_ch_rel, max_ch_abs, n, &(L_b), zt1_b, zt2_b, &UL, Uzt1, Uzt2); 
            fprintf(stderr, " ch_rel = %12.8f ch_abs = %12.8f\n", ch_rel, ch_abs);
            stm_compute_mismatch(FALSE, FALSE, n, p1, p2, &(L_b), zt1_b, zt2_b, &S_b, NULL, NULL, NULL);
            tries++; 
            if (S_b >= S) 
              { double a = (S_b - S + h*b)/(h*h);
                double h_min_est = b/a/2.0;
                if (h_min_est >= 0.9375*h) 
                  { h = 0.9375*h; fprintf(stderr, ">"); }
                else if (h_min_est <= 0.0625*h) 
                  { h = 0.0625*h; fprintf(stderr, "<"); }
                else
                  { h = h_min_est; fprintf(stderr, "|"); }
              }
            else if (tries > max_tries) 
              { fprintf(stderr, "*** unimin: %d attempts - giving up ***\n", tries);
                exit(1);
              }
          }
        fprintf(stderr, "×");
        /*  Move to improved state: */
        (*L_p) = L_b;
        for (k = 0; k < n; k++) { zt1[k] = zt1_b[k]; zt2[k] = zt2_b[k]; }
        /*  S_a = S_b; */
        /*  DS_DL_a = DS_DL_b; */
        /*  t = DS_Dzt1_a; DS_Dzt1_a = DS_Dzt1_b; DS_Dzt1_b = t; */
        /*  t = DS_Dzt2_a; DS_Dzt2_a = DS_Dzt2_b; DS_Dzt2_b = t; */
      }
    }
}

void select_direction(
  int n, r2_t *p1, r2_t *p2, 
  stereo_params_t *L_p, double zt1[], double zt2[],
  stereo_params_t *DS_DL_p, double DS_Dzt1[], double DS_Dzt2[],
  stereo_params_t *UL_p, double Uzt1[], double Uzt2[])
  {
    double D1, D2 = 0;
    double x1, y1, z1, x2, y2, z2;
    int i,j,k;
    /*  Compute the approximate mean spread of the points on each side: */
    for (k = 0; k < n; k++)
      { x1 = p1[k].c[0]; y1 = p1[k].c[1]; z1 = zt1[k];
        D1 += x1*x1 + y1*y1 + z1*z1;
        x2 = p2[k].c[0]; y2 = p2[k].c[1]; z2 = zt2[k];
        D2 += x2*x2 + y2*y2 + z2*z2;
      }
    D1 = sqrt(D1/n); D2 = sqrt(D2/n);
    /*  Use the gradient direction, but deflate {R} and {m} by {D1,D2}: */
    { double m = DS_DL_p->m;
      double fR = (D1*m + D2/m);
      double fm = (D1 + D2/(m*m));
      (*UL_p) = (*DS_DL_p);
      for (i = 0; i < 3; i++)
        { for (j = 0; j < 3; j++)
            { UL_p->R.c[i][j] /= fR; }
        }
      UL_p->m /= fm;
      for (k = 0; k < n; k++)
        { Uzt1[k] = DS_Dzt1[k]; Uzt2[k] = DS_Dzt2[k]; }
    }
  }

