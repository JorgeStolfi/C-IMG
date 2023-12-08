// stereo-mismatch.h - Procedures to compute the mismatch between two stereo images
// Last edited on 2015-06-12 21:21:38 by stolfilocal

#ifndef stereo_mismatch_H
#define stereo_mismatch_H

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <r2.h>
#include <r3.h>
#include <r3x3.h>

#include <stereo_basics.h>

typedef struct stereo_params_t {
  double d;      // Distance from camera to image plane, in pixels.
  r3x3_t R;      // Rotation matrix.
  double m;      // Scale factor from image 1 to image 2.
  r3_t v;        // Displacement (after scaling).
} stereo_params_t;
  /* Parameters of the stereo correspondence. The correspondence is
    defined by the product of 4×4 homogeneous matrices
      |
      | {Q = P^{-1} R S T P}
      | 
    where the matrices are derived from the parameters {L=(d,R,m,v)} thusly:
      |
      | {P} = conical perspective matrix = {[[d 0 0 0] [0 d 0 0] [0 0 d 0] [-1 0 0 d]]}.
      | {R} = a rotation matrix, defined by a 3×3 orthonormal matrix {R}.
      | {S} = a uniform scaling matrix, with scale factor {m}.
      | {T} = a translation by the 3-vector {v}.
      |
    The matrix {Q} maps the homogeneous Perspective coordinates
    {[wp1,xp1,yp1,zp1]} of a point in one image to the Perspective
    coordinates {[wp2,xp2,yp2,zp2]} of the same point in the second image.

    In the Perspective Coordinate System (PCS), the coordinates of a
    scene point are those resulting from the conical transformation that
    sends the camera from {(0,0,d)} to {(0,0,+oo)}, before the
    orthogonal projection onto the image. In particular, {(xp/wp,yp/wp)}
    are the Cartesian coordinates of the point's projection onto the
    image, relative to the image center.

    More generally, the PCS coordinates are related to the true
    Cartesian coordinates {(Xt,Yt,Zt)} in the Image Coordinate System
    (ICS) by the formula {(xp/wp,yp/wp,zp/wp) = f * (Xt,Yt,Zt)},
    where {f} is the perspective scaling factor, {f = d/(d - Zt)}.
    The camera has height {Zp = +oo} in the PCS, {Zt = d} in the ICS.

    This formulation assumes that in both coordinate systems, PCS and
    ICS, the origin is the perpendicular projection of the camera onto
    the image plane (usually, the center of the image). It also assumes
    that the images are isotropic, distortion-free, and were obtained
    with the same camera lens and and zoom setting (i.e. with the same
    {d} parameter). */

void stm_guess_translation(int n, r2_t *p1, r2_t *p2, 
  stereo_params_t *L_p, double zt1[], double zt2[]);
  /* Computes an initial guess for the translation term {L_p->v}, given 
    the apparent (PCS) coordinates of the corresponding pairs of points 
    {p1[k] <-> p2[k]}, their guessed ICS heights {zt1[k],zt2[k]}, and
    the other correspondence parameters in {*L_p}. */

void stm_invert_XYZ_perspective(r3_t *p_p, stereo_params_t *L_p, r3_t *t_p);
  /* Computes the true ICS Cartesian coordinates {*t_p = (Xt,Yt,Zt)}
    of the point whose PCS Cartesian coordinates are {*p_p = (Xp,Yp,Zp)}. */

void stm_invert_XY_perspective(r2_t *p_p, double zt, stereo_params_t *L_p, r2_t *t_p);
  /* Given the apparent (PCS) Cartesian coordinates {*p_p = (xp,yp)} of a point,
    and its true height {zt} from the image plane, computes its true image-relative
    (ICS) Cartesian coordinates {*t_p = (xt,yt)}. */

void stm_compute_mismatch
  ( bool_t grad, 
    bool_t verbose, 
    int n, r2_t *p1, r2_t *p2, 
    stereo_params_t *L_p, double zt1[], double zt2[],
    double *S_p, 
    stereo_params_t *DS_DL_p, double DS_Dzt1[], double DS_Dzt2[]
  );
  /* Computes the total quadratic mismatch {*S_p} between the predicted
    and actual coordinates of the paired points. If {grad = true}, also
    computes the derivatives of {*S_p} with respect to all
    correspondence parameters in {*L_p} and to the presumed point
    heights {zt1[k]}, {zt2[k]}.

    The mismatch is defined as the sum of the squared differences
    between the PCS coordinates {(xp,yp)}, because they are observed
    quantites, and between the ICS heights {zt}, because they are more
    stable than the PCS heights {zp}. The coorrespondence transformation
    is applied both directions (from image 1 to image 2, and vice-versa)
    and the two mismatches are added. */
  
double stm_gradient_dot(int n, 
  stereo_params_t *UL_p, double Uzt1[], double Uzt2[], 
  stereo_params_t *VL_p, double Vzt1[], double Vzt2[]);
  /* Computes the dot product of two vectors in parameter space. */

double stm_gradient_norm_sqr(int n, 
  stereo_params_t *UL_p, double Uzt1[], double Uzt2[]);
  /* Computes the square of the norm of a vector in parameter space. */

void stm_modify_params(double h, 
  double ch_rel_p[], double ch_abs_p[], 
  double max_ch_rel, double max_ch_abs,
  int n, 
  stereo_params_t *L_p, double zt1[], double zt2[],
  stereo_params_t *UL_p, double Uzt1[], double Uzt2[]);
  /* Changes the parameters {(*L_p, *zt1, *zt2)} by (approximately) {h}
    times the negative of the direction {(*UL_p, *Uzt1, *Uzt2)},
    ensuring that the {d} and {m} components remain positive, and {R}
    remains orthonormal.

    Also stores in {*ch_rel_p} the relative change (log ratio) in the
    parameters {d}, {m} and {R}. Stores in {*ch_abs_p} the absolute
    change (difference) in the parameters {v}, {zt1}, {zt2}. The maximum
    change in each component is a factor of {exp(max_ch_rel)} for {d},
    {m} and {R}, {max_ch_abs} for {v}, {zt1}, {zt2}. */

void stm_print_params(FILE *f, stereo_params_t *L_p);
  /* Prints the correspondence parameters {*L_p} to the file *f. */

void stm_print_corresp(FILE *f, int n, stereo_params_t *L_p, double zt1[], double zt2[]);
  /* Prints the correspondence parameters {*L_p} and the estimated
    heights {zt1,zt2} to the file *f. */

void stm_test_gradient(int n, r2_t *p1, r2_t *p2, stereo_params_t *L_p, double zt1[], double zt2[]);
  /* Performs a test of the gradient code in {stm_compute_mismatch} 
    by comparing its result to that of numerical derivation. */

#endif
