#define PROG_NAME "image_stitch"
#define PROG_DESC "Finds projective map between two images, given corresponding points"
#define PROG_VERS "1.0"

// Last edited on 2024-12-21 14:00:51 by stolfi

#define image_stitch_C_COPYRIGHT \
    "© 2002 by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    -image1 {LX1} {LY1} {HX1} {HY1} \\\n" \
  "    -image2 {LX2} {LY2} {HX2} {HY2} \\\n" \
  "    [ -maxIter {MAXITER} ] \\\n" \
  "    [ -maxErr {MAXERR} ] \\\n" \
  "    [ -verbose ] \\\n" \
  "    -outPrefix {OPREF} \\\n" \
  "    < {PAIRFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "   Reads a file of corresponding points in two images.  Computes" \
  " a projective map for each image that produce a good match for" \
  " those points and minimizes the distortion.\n" \
  "\n" \
  "   The two projective maps will be inverses of each" \
  " other, so that the merged image will be a meet-in-the-middle" \
  " compromise.\n" \
  "\n" \
  "   The input file must have one point pair per" \
  " line, in the format \"( {X1} {Y1} ) = ( {X2} {Y2} )\" where" \
  " {(X1,Y1)} is a point in the domain of image 1, and {(X2,Y2)} is" \
  " a point in the domain of image 2.  The coordinates may be expressed" \
  " in any units and any coordinate system; the matrices will be" \
  " expressed in the same units and system.\n" \
  "\n" \
  "OUTPUT FILES\n" \
  "  The program writes to \"{OPREF}-{K}-matrix.txt\" the matrix to be used for" \
  " image number {K} (\"1\" or \"2\").\n" \
  "\n" \
  "  It also writes to \"{OPREF}-pairs.txt\" the original and mapped point" \
  " pairs.  The file will have one line for each point pair, with eight" \
  " fields: old X, old Y, new X and new Y for the point in image 1, and" \
  " then the same for the point in image 2.\n" \
  "\n" \
  "  It also writes to \"{OPREF}-{K}-outline.txt\" the outlines of image {K} mapped" \
  " by the corresponding direct matrix.  The file has five lines, with" \
  " the X and Y coordinates of each corner in each line, in counterclockwise" \
  " order around the original domain, with the fifth point being a copy of" \
  " the first.\n" \
  "\n" \
  "  If non-linear optimization was used, also writes to \"{OPREF}-f2-plot.txt\" a" \
  " plot of the goal function (mean quadratic point mismatch) in the" \
  " neighborhood of the computed optimum matrices.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -image1 {LX1} {LY1} {HX1} {HY1}\n" \
  "  -image2 {LX2} {LY2} {HX2} {HY2}\n" \
  "    These mandatory arguments specify the nominal domains of the" \
  " two images, namely {[LXi _ HXi]×[LYi _ HYi]}, in" \
  " the same units are the point coordinates.  They are used" \
  " to compute the amount of distortion; the data points need not be" \
  " inside these nominal domains.  Typically, {LXi} and {LYi} are" \
  " zero, {HXi} is the width of image {i}, and {HYi} is its height.\n" \
  "\n" \
  "  -maxIter {MAXITER}\n" \
  "    This optional argument is the number of iterations for" \
  " the non-linear optimizitation.  The default is \"-maxIter 5\".\n" \
  "\n" \
  "  -maxErr {MAXERR}\n" \
  "    This optional argument is the convergence criterion for non-linear" \
  " optimization.  The optimization will stop when the roo mean square" \
  " mismatch of the mapped points is less than {MAXERR}.  If {MAXERR}" \
  " is zero (the default), the" \
  " optimization will continue for {MAXITER} iterations.\n" \
  "\n" \
  "  -outPrefix {OPREF}\n" \
  "    This mandatory argument specifies the common prefix of all output files.\n" \
  "\n" \
  "  -verbose\n" \
  "    If present, this option requests debugging printouts to {stderr}.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  Hear also.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created nov/2011 by J.Stolfi.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " image_stitch_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <r2x2.h>
#include <r3.h>
#include <r3x3.h>
#include <rn.h>
#include <hr2.h>
#include <r2_aff_map.h>
#include <jsfile.h>
#include <sve_minn.h>
#include <minn_plot.h>
#include <interval.h>
#include <interval_io.h>
#include <box.h>
#include <argparser.h>

typedef struct options_t 
  { r2_t L1, H1;          /* Nominal domain corners of first image. */
    r2_t L2, H2;          /* Nominal domain corners of second image. */
    char *outPrefix;      /* Prefix for output file names. */
    bool_t verbose;       /* TRUE to print debugging info. */
    /* Parameters of the non-lienar optimizer: */
    int32_t maxIter;      /* Max number of iterations. */
    double maxErr;        /* Convergence criterion. */
  } options_t;

options_t *image_stitch_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */

void image_stitch_read_data(r2_vec_t *p1, r2_vec_t *p2, bool_t verbose);
/* Reads from {stdin} a list of point pairs {p1[i],p2[i]},
  where {p1[i]} is in image 1 and {p2[i]} is in image 2.
  Each pair should be in the format "( h1 v1 ) = ( h2 v2 )" where 
  {h1,h2} are column indices and {v1,v2} are row indices,
  both counted from 0.  The vectors {p1,p2} are allocated
  by the procedure. */

hr2_pmap_t image_stitch_initial_pmap_guess(r2_vec_t *p1, r2_vec_t *p2, bool_t verbose);
  /* Computes an affine matrix {M} such that {p1} mapped TWICE by {M} 
    matches {p2}, approximately. */

hr2_pmap_t image_stitch_optimize_pmap
  ( hr2_pmap_t *M0, 
    r2_t *L1, 
    r2_t *H1, 
    r2_vec_t *p1, 
    r2_t *L2, 
    r2_t *H2, 
    r2_vec_t *p2, 
    int32_t maxIter,
    double maxErr,
    char *outPrefix, 
    bool_t verbose
  );
  /* Computes the projective matrix {M} such that {p1} mapped TWICE by {M} 
    is as close as possible to {p2}, by linear optimization. 
    Also writes plottable samples of the quadratic goal function in 
    the neighborhood of the computed optimum to 
    file "{outPrefix}-f2-plot.txt". */

void image_stitch_map_point(r2_t *p, r3x3_t *M, r2_t *q);
  /* Maps the Cartesian point {p} through the projective map with 
    homogeneous matrix {M}, stores result in {*q}. Returns {(INF,INF)} 
    if the resulting point is at infinity or beyond. */

void image_stitch_bar(r2_vec_t *p, r2_t *bar);
  /* Computes the barycenter {bar} of the point set {*p}. */

double image_stitch_mean_err_sqr(r2_vec_t *p1, r3x3_t *M1, r2_vec_t *p2, r3x3_t *M2);
  /* Cmputes the mean squared error between the points of {p1} mapped by 
    {M1} and the corresponding points of {p2} maped by {M2}. */

double image_stitch_deform_sqr(r2_t *Lp, r2_t *Hp, r3x3_t *M);
  /* Computes a quadratic deformation term for the axis-aligned rectangle with  
    corners {Lp} and {Hp} when mapped through the projective map {M}. 
    The result is 0 for Euclidean maps, and positive for other maps. */

void image_stitch_show_pmap(hr2_pmap_t *M, char *tag);
  /* Writes (to stderr) the direct and inverse matrices of {M}. */

void image_stitch_show_pmap_and_square(hr2_pmap_t *M, char *tag);
  /* Writes (to stderr) the direct and inverse matrices of {M} and of its 
    square {M*M}. */

void image_stitch_write_matrix(char *outPrefix, int32_t K, r3x3_t *M);
  /* Writes the matrix {M} for image {K} (1 or 2)
    as file "{outPrefix}-{K}-matrix.txt". */

void image_stitch_compute_mapped_outline(r2_t *L, r2_t *H, r3x3_t *M, int32_t nc, r2_t C[]);
  /*Stores in {C[0..3]} the corners of an image, mapped by {M}. Assumes {L} and {H} are the inferior
    and superior corners of the unmapped image domain. */

void image_stitch_write_points(char *outPrefix, int32_t K, int32_t nc, r2_t C[]);
  /* Writes to file "{outPrefix}-{K}-outline.txt" the points {C[0..nc-1]}, assumed
    to be the corners of a polygon.  The point {C[0]} is repeated at the end
    to close the polygon. */

void image_stitch_check_matrices(char *outPrefix, r2_vec_t *p1, r3x3_t *M1, r2_vec_t *p2, r3x3_t *M2);
  /* Compares (to stderr) the points of {p1} mapped by 
    {M1} to the corresponding points of {p2} maped by {M2}.
    If {outPrefix} is NULL writes the report to {stderr},
    else writes to file "{outPrefix}-pairs.txt". */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = image_stitch_parse_options(argc, argv);
    
    /* Read the point pairs: */
    r2_vec_t p1;    /* Reference points in image 1. */
    r2_vec_t p2;    /* Reference points in image 2. */
    image_stitch_read_data(&p1, &p2, o->verbose);
    int32_t np = p1.ne;
    assert(np == p2.ne);
    if (np < 3) { fprintf(stderr, "too few point pairs (%d) for affine mapping\n", np); exit(1); }
    
    hr2_pmap_t M0 = image_stitch_initial_pmap_guess(&p1, &p2, o->verbose);
                                        
    if (o->verbose)
      { /* Checking: */
        image_stitch_show_pmap_and_square(&M0, "affine");
        fprintf(stderr, "affine mapped points:\n");
        image_stitch_check_matrices("", &p1, &(M0.dir), &p2, &(M0.inv));
      }

    /* Nonlinear optimization: */
    hr2_pmap_t M = image_stitch_optimize_pmap
      ( &M0, 
        &(o->L1), &(o->H1), &p1, 
        &(o->L2), &(o->H2), &p2, 
        o->maxIter, o->maxErr, 
        o->outPrefix, o->verbose
      );

    /* Check and show: */
    image_stitch_show_pmap_and_square(&M, "final");
    fprintf(stderr, "final mapped points:\n");
    image_stitch_check_matrices("", &p1, &(M.dir), &p2, &(M.inv));
    
    /* Compute mapped outlines: */
    int32_t nc = 4;
    r2_t C1[4], C2[4]; 
    image_stitch_compute_mapped_outline(&(o->L1), &(o->H1), &(M.dir), nc, C1);
    image_stitch_compute_mapped_outline(&(o->L2), &(o->H2), &(M.inv), nc, C2);
    
    /* Output and checking: */
    image_stitch_write_matrix(o->outPrefix, 1, &(M.dir));
    image_stitch_write_matrix(o->outPrefix, 2, &(M.inv));
    image_stitch_write_points(o->outPrefix, 1, nc, C1);
    image_stitch_write_points(o->outPrefix, 2, nc, C2);
    image_stitch_check_matrices(o->outPrefix, &p1, &(M.dir), &p2, &(M.inv));
    
    return 0;
  }

void image_stitch_read_data(r2_vec_t *p1, r2_vec_t *p2, bool_t verbose)
  {
    if (verbose) { fprintf(stderr, "--- enter image_stitch_read_data ---\n"); }
    (*p1) = r2_vec_new(20);
    (*p2) = r2_vec_new(20);
    int32_t np = 0; /* Number of data pairs. */

    while (TRUE)
      { double x1, y1, x2, y2;
        int32_t nscan = fscanf(stdin, " ( %lf %lf ) = ( %lf %lf )", &x1, &y1, &x2, &y2);
        if (nscan <= 0) { break; }
        if (nscan != 4) { fprintf(stderr, "line %d: bad input format (nscan = %d)\n", np+1, nscan); exit(1); }
        fprintf(stderr, " ( %6.1f %6.1f ) = ( %6.1f %6.1f )\n", x1,y1,x2,y2);
        r2_vec_expand(p1, np); p1->e[np] = (r2_t){{x1, y1}}; 
        r2_vec_expand(p2, np); p2->e[np] = (r2_t){{x2, y2}}; 
        np++;
      }
    r2_vec_trim(p1, np); 
    r2_vec_trim(p2, np);
    if (verbose) { fprintf(stderr, "--- exit image_stitch_read_data ---\n"); }
  }

hr2_pmap_t image_stitch_initial_pmap_guess(r2_vec_t *p1, r2_vec_t *p2, bool_t verbose)
  {
    int32_t np = p1->ne;
    assert(np == p2->ne);
    
    if (verbose) { fprintf(stderr, "--- enter image_stitch_initial_pmap_guess ---\n"); }
    r2_aff_map_t A = r2_aff_map_from_point_pairs(np, p1->e, p2->e, NULL);
    hr2_pmap_t M = hr2_pmap_from_mat_and_disp(&(A.mat), &(A.disp));
    
    /* Split matrix into translation vector {v} and linear submatrix {S}: */
    assert(M.dir.c[0][0] == 1);
    assert(M.dir.c[1][0] == 0);
    assert(M.dir.c[2][0] == 0);
    r2_t v = (r2_t){{ M.dir.c[0][1], M.dir.c[0][2] }};
    r2x2_t S;
    S.c[0][0] = M.dir.c[1][1];
    S.c[0][1] = M.dir.c[1][2];
    S.c[1][0] = M.dir.c[2][1];
    S.c[1][1] = M.dir.c[2][2];
    
    /* Obtain square root {R} of {S}: */
    r2x2_t R; r2x2_sqrt(&S, &R);

    /* Compute the translation vetor {u} so that {(p*R +u)R + u = p*R^2 + v}: */
    r2x2_t N = R;
    N.c[0][0] += 1;
    N.c[1][1] += 1;
    r2x2_inv(&N, &N);
    r2_t u; 
    r2x2_map_row(&v, &N, &u);
    
    /* Pack {R} and {u} as a 3x3 homogeneous affine matrix {M}: */
    hr2_pmap_t Q;
    Q.dir.c[0][0] = 1.0;
    Q.dir.c[0][1] = u.c[0];
    Q.dir.c[0][2] = u.c[1];

    Q.dir.c[1][0] = 0.0;
    Q.dir.c[1][1] = R.c[0][0];
    Q.dir.c[1][2] = R.c[0][1];

    Q.dir.c[2][0] = 0.0;
    Q.dir.c[2][1] = R.c[1][0];
    Q.dir.c[2][2] = R.c[1][1];
    
    r3x3_inv(&(Q.dir), &(Q.inv));
    if (verbose) { fprintf(stderr, "--- exit image_stitch_initial_pmap_guess ---\n"); }
    return Q;
  }
        

hr2_pmap_t image_stitch_optimize_pmap
  ( hr2_pmap_t *M0, 
    r2_t *L1,
    r2_t *H1,
    r2_vec_t *p1, 
    r2_t *L2,
    r2_t *H2,
    r2_vec_t *p2, 
    int32_t maxIter,
    double maxErr,
    char *outPrefix,
    bool_t verbose
  )
  {
    bool_t debug = TRUE;

    int32_t np = p1->ne;
    assert(np == p2->ne);
    
    if (verbose) { fprintf(stderr, "--- enter image_stitch_optimize_pmap ---\n"); }
    if (np < 4)
      { fprintf(stderr, "too few point pairs (%d) for projective map, using affine map\n", np); 
        return (*M0);
      }
    
    int32_t nx = 8; /* Number of packed parameters. */

    auto hr2_pmap_t compute_S_map(void);
      /* Computes the {S} map used by {unpack_parameters}. */
         
    hr2_pmap_t S = compute_S_map();  /* Pmap from image1's domain to canonical square. */

    auto void pack_parameters(hr2_pmap_t *M, int32_t nx, double x[]);
      /* Packs the homogeneous projective map {M} to the parameter
        vector {x[0..nx-1]}, as the images of the four corners 
        of image 1. Requires {nx==8}. */
         
    auto void unpack_parameters(int32_t nx, double x[], hr2_pmap_t *M);
      /* Unpacks the parameter vector {x[0..nx-1]} to the 
        homogeneous map matrix {M}.  The matrix
        will have {M[0,0] == 1}. Requires {nx==8}. */
         
    auto double goalf(int32_t nx, double x[]);
      /* Computes the mean squared distance between the positions of the
         mapped points {p1*M.dir} and {p2*M.inv}, given the packed
         parameters {x[0..nx-1]}. */
    
    auto bool_t is_ok(int32_t iter, int32_t nx, double x[], double Fx, double dist, double step, double radius);
      /* Returns true if the squared error {Fx} is small enough. */
    
    double esq = image_stitch_mean_err_sqr(p1, &(M0->dir), p2, &(M0->inv));
    double dMax = 20*sqrt(esq);      /* Search radius around initial guess. */
    bool_t dBox = FALSE;             /* Search in ball, not in box. ??? */
    double rIni = 0.250*dMax;        /* Initial probe radius. */
    double rMin = 0.5;               /* Minimum probe radius. */
    double rMax = 0.500*dMax;        /* Maximum probe radius. */
    double minStep = 0.1*maxErr;        /* Stopping criterion. */
    sign_t dir = -1;                 /* Look for minimum. */

    if (verbose) 
      { fprintf(stderr, "initial mean error squared = %22.16e\n", esq);
        fprintf(stderr, "estimated distance from optimum = %13.6f\n", dMax);
        fprintf(stderr, "probe radius = %13.6f [ %13.6f _ %13.6f ]\n", rIni, rMin, rMax);
      }
      
    double x[nx]; /* Initial guess and final optimum parameters. */ 
    pack_parameters(M0, nx, x);
    if (verbose) 
      { char *fname = NULL;
        char *fname = jsprintf("%s-f2-plot.txt", outPrefix);
        FILE *wr = open_write(fname, TRUE);
        for (uint32_t ku = 0;  ku < nu; ku++)
          { double *uk = &(U[ku*nz]);
            minn_plot_1D_gnuplot(wr, nx, goalf, x, 20, rIni);
          }
        fclose(wr);
        free(fname);
      }
    if (verbose) { fprintf(stderr, "optimizing\n"); }
    double Fx = goalf(nx, x);
    if (verbose) { fprintf(stderr, "initial rms error = %13.6f\n", Fx); }
    double ctr[nx]; rn_copy(nx, x, ctr);
    bool_t sve_debug = verbose;
    bool_t sve_debug_probes = FALSE; 
    sve_minn_iterate
      ( nx, &goalf, &is_ok, NULL, 
        x, &Fx, dir, 
        ctr, dMax, dBox, rIni, rMin, rMax, 
        minStep, maxIter, sve_debug, sve_debug_probes
      );
    if (verbose) { fprintf(stderr, "final rms error = %13.6f\n", Fx); }
    
    /* Unpack projective map: */
    hr2_pmap_t M; /* Optimized map. */
    unpack_parameters(nx, x, &M);
    if (verbose) { fprintf(stderr, "--- exit image_stitch_optimize_pmap ---\n"); }
    return M;
 
    /* --- internal procs ----------------------------------------------------- */

    hr2_pmap_t compute_S_map(void)
      { hr2_point_t hp[4];
        int32_t ix, iy;
        for (ix = 0; ix < 2; ix++)
          { double px = (ix == 0 ? L1->c[0] : H1->c[0]);
            for (iy = 0; iy < 2; iy++)
              { double py = (iy == 0 ? L1->c[1] : H1->c[1]);
                hp[2*iy+ix] = (hr2_point_t){{{ 1, px, py }}};
              }
          }
        hr2_pmap_t S = hr2_pmap_from_four_points(&(hp[0]), &(hp[1]), &(hp[2]), &(hp[3]));
        return hr2_pmap_inv(&S);
      }
    
    void pack_parameters(hr2_pmap_t *M, int32_t nx, double x[])
      { assert(nx == 8);
        int32_t ix, iy;
        int32_t k = 0;
        for (ix = 0; ix < 2; ix++)
          { double px = (ix == 0 ? L1->c[0] : H1->c[0]);
            for (iy = 0; iy < 2; iy++)
              { double py = (iy == 0 ? L1->c[1] : H1->c[1]);
                r2_t p = (r2_t){{ px, py }};
                r2_t q; image_stitch_map_point(&p, &(M->dir), &q);
                x[k] = q.c[0]; k++;
                x[k] = q.c[1]; k++;
              }
          }
      }
    
    void unpack_parameters(int32_t nx, double x[], hr2_pmap_t *M)
      { assert(nx == 8);
        int32_t ix, iy;
        int32_t k = 0;
        hr2_point_t hp[4];
        for (ix = 0; ix < 2; ix++)
          { for (iy = 0; iy < 2; iy++)
              { double qx = x[k]; k++;
                double qy = x[k]; k++;
                hp[2*iy + ix] = (hr2_point_t){{{ 1, qx, qy }}};
              }
          }
        hr2_pmap_t R = hr2_pmap_from_four_points(&(hp[0]), &(hp[1]), &(hp[2]), &(hp[3]));
        (*M) = hr2_pmap_compose(&S, &R); 
      }
    
    double goalf(int32_t nx, double x[])
      {
        if (sve_debug) { rn_gen_print(stderr, nx, x, "%8.1f", "\n    [ ", " ", " ]\n"); }
        hr2_pmap_t M;
        unpack_parameters(nx, x, &M);
        /* Mean square error on point lists: */
        double esq = hr2_pmap_mismatch_sqr(&M, np, p1, p2, w);
        /* Deformation penalty: */
        double dsq1 = image_stitch_deform_sqr(L1, H1, &(M.dir));
        double dsq2 = image_stitch_deform_sqr(L2, H2, &(M.inv));
        double F = esq + 0.0001*(dsq1 + dsq2);
        if (sve_debug) 
          { fprintf(stderr, "    mean squared error = %13.6e", esq); 
            fprintf(stderr, " squared deform = %13.6e %13.6e", dsq1, dsq2);
            fprintf(stderr, " function = %13.6e\n", F);
          }
        return F;
        
      }
      
    void report(hr2_pmap_t *M, double F)
      { 
      
      }
      
    bool_t is_ok(int32_t iter, int32_t nx, double x[], double Fx, double dist, double step, double radius)
      {
        return Fx < maxErr*maxErr;
      }
    
  }

void image_stitch_map_point(r2_t *p, r3x3_t *M, r2_t *q)
  {
    r3_t hp = (r3_t){{ 1, p->c[0], p->c[1] }};
    r3_t hq;
    r3x3_map_row(&hp, M, &hq);
    double w = hq.c[0];
    double m = fmax(fabs(hq.c[1]), fabs(hq.c[2]));
    if (w <= m*1e-200) 
      { (*q) =  (r2_t){{ INF, INF }}; }
    else
      { (*q) =  (r2_t){{ hq.c[1]/w, hq.c[2]/w }}; } 
  }

void image_stitch_show_pmap_and_square(hr2_pmap_t *M, char *tag)
  {
    image_stitch_show_pmap(M, tag);
    char *tagsqr = NULL;
    char *tagsqr = jsprintf("%s squared", tag);
    hr2_pmap_t M2 = hr2_pmap_compose(M,M);
    image_stitch_show_pmap(&M2, tagsqr);
    free(tagsqr);
  }

void image_stitch_show_pmap(hr2_pmap_t *M, char *tag)
  {
    fprintf(stderr, "  %s matrix (direct):\n", tag);
    r3x3_gen_print(stderr, &(M->dir), "%13.6e", "", "\n", "\n", "    [ ", " ", " ]");

    fprintf(stderr, "  %s matrix (inverse):\n", tag);
    r3x3_gen_print(stderr, &(M->inv), "%13.6e", "", "\n", "\n", "    [ ", " ", " ]");
  }
    
void image_stitch_check_matrices(char *outPrefix, r2_vec_t *p1, r3x3_t *M1, r2_vec_t *p2, r3x3_t *M2)
  {
    int32_t np = p1->ne;
    assert(np == p2->ne);
    
    bool_t delims = TRUE;  /* Shall we print delimiters and arrows? */
    
    /* Choose the output file: */
    FILE *wr = stderr;     
    if ((outPrefix != NULL) && (strlen(outPrefix) != 0))
      { char *fname = NULL;
        char *fname = jsprintf("%s-pairs.txt", outPrefix);
        wr = open_write(fname, TRUE);
        delims = FALSE;
        free(fname);
      }
    
    if (delims) 
      { fprintf(wr, "columns: ( X1 Y1 ) -> ( X1' Y1' )  ( XERR YERR ) ERR  ( X2' Y2' ) <- ( X2 Y2 )\n"); }  
    else
      { fprintf(wr, "# columns: X1 Y1  X1' Y1'  X2 Y2  X2' Y2'  XERR YERR  ERR\n"); }  
    
    int32_t k;
    double sum2 = 0;
    double emax = -INF;
    for(k = 0; k < np; k++)
      { r2_t *p1k = &(p1->e[k]);
        r2_t q1k; image_stitch_map_point(p1k, M1, &q1k);
        r2_t *p2k = &(p2->e[k]);
        r2_t q2k; image_stitch_map_point(p2k, M2, &q2k);
        r2_t d; r2_sub(&q1k, &q2k, &d);
        double e2 = r2_dist_sqr(&q1k, &q2k);
        double e = sqrt(e2);
                  
        if (delims)
          { 
            r2_gen_print(wr, p1k, "%8.2f", "( ", " ", " )");
            fprintf(wr, " -> ");
            r2_gen_print(wr, &q1k, "%8.2f", "( ", " ", " )");
            
            fprintf(wr, "  ");
            
            r2_gen_print(wr, &d, "%+6.2f", "( ", " ", " )");
            fprintf(wr, " ");
            fprintf(wr, "%6.2f", e);
            
            fprintf(wr, "  ");
            
            r2_gen_print(wr, &q2k, "%8.2f", "( ", " ", " )");
            fprintf(wr, " <- ");
            r2_gen_print(wr, p2k, "%8.2f", "( ", " ", " )");
          }
        else
          {
            r2_gen_print(wr, p1k, "%12.6e", "", " ", "");
            fprintf(wr, "  ");
            r2_gen_print(wr, &q1k, "%12.6e", "", " ", "");
            
            fprintf(wr, "  ");
            
            r2_gen_print(wr, p2k, "%12.6e", "", " ", "");
            fprintf(wr, "  ");
            r2_gen_print(wr, &q2k, "%12.6e", "", " ", "");
            
            fprintf(wr, "  ");
            
            r2_gen_print(wr, &d, "%+12.6e", "", " ", "");
            fprintf(wr, " ");
            fprintf(wr, "%12.6e", e);
          }
          
        fprintf(wr, "\n");
        
        sum2 += e2;
        emax = fmax(emax, e);
      }
    double erms = sqrt(sum2/np);
    if (delims)
      { fprintf(wr, "root mean squared error = %8.4f\n", erms);
        fprintf(wr, "max error               = %8.4f\n", emax);
      }
    if (wr != stderr) { fclose(wr); }
  }

void image_stitch_write_matrix(char *outPrefix, int32_t K, r3x3_t *M)
  {
    assert((K == 1) || (K == 2));
    char *fname = jsprintf("%s-%d-matrix.txt", outPrefix, K);
    FILE *wr = open_write(fname, TRUE);
    r3x3_gen_print(wr, M, "%24.15e", "", "", "", "", " ", "\n");
    fclose(wr);
    free(fname);
  }

void image_stitch_compute_mapped_outline(r2_t *L, r2_t *H, r3x3_t *M, int32_t nc, r2_t C[])
  {
    assert(nc == 4);

    int32_t kc = 0; /* Number of corners computed. */
    
    auto void comp_corner(double X, double Y);
      /* Maps the point {(X,Y)} by the projective map with matrix {M} and 
        appends the result to {C[nc]}, inclrementing {nc}. */
    
    comp_corner(L->c[0], L->c[1]); 
    comp_corner(H->c[0], L->c[1]); 
    comp_corner(H->c[0], H->c[1]); 
    comp_corner(L->c[0], H->c[1]); 
    
    assert(kc == nc);
    
    return;
    
    void comp_corner(double X, double Y)
      { r2_t p = (r2_t){{ X, Y }};
        r2_t q; image_stitch_map_point(&p, M, &q);
        assert(kc < nc);
        C[kc] = q;
        kc++;
      }
  }

void image_stitch_write_points(char *outPrefix, int32_t K, int32_t nc, r2_t C[])
  {
    assert((K == 1) || (K == 2));
    
    char *fname = jsprintf("%s-%d-outline.txt", outPrefix, K);
    FILE *wr = open_write(fname, TRUE);
    for (uint32_t kc = 0;  kc <= nc; kc++)
      { r2_t *q = &(C[kc % nc]);
        r2_gen_print(wr, q, "%12.6e", "", " ", "\n");
      }
    fclose(wr);
    return;
  }

options_t *image_stitch_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    /* Parse keyword-based arguments: */
    argparser_get_keyword(pp, "-image1");
    o->L1.c[0] = argparser_get_next_double(pp, -100000, +100000);
    o->L1.c[1] = argparser_get_next_double(pp, -100000, +100000);
    o->H1.c[0] = argparser_get_next_double(pp, -100000, +100000);
    o->H1.c[1] = argparser_get_next_double(pp, -100000, +100000);
                                             
    argparser_get_keyword(pp, "-image2");    
    o->L2.c[0] = argparser_get_next_double(pp, -100000, +100000);
    o->L2.c[1] = argparser_get_next_double(pp, -100000, +100000);
    o->H2.c[0] = argparser_get_next_double(pp, -100000, +100000);
    o->H2.c[1] = argparser_get_next_double(pp, -100000, +100000);
    
    if (argparser_keyword_present(pp, "-maxIter"))
      { o->maxIter = (int32_t)argparser_get_next_int(pp, 0, 100000); }
    else
      { o->maxIter = 5; }

    if (argparser_keyword_present(pp, "-maxErr"))
      { o->maxErr = argparser_get_next_double(pp, 0, +INF); }
    else
      { o->maxErr = 0.0; }
    
    if (argparser_keyword_present(pp, "-writeStitched"))
      { o->writeStitched = argparser_get_next_bool(pp); }
    else
      { o->writeStitched = FALSE; }

    argparser_get_keyword(pp, "-outPrefix");    
    o->outPrefix = argparser_get_next_non_keyword(pp);
    
    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }


