#define PROG_NAME "pnmalign"
#define PROG_DESC "find optimum relative alignment of pbm/ppm/pgm files"
#define PROG_VERS "1.0"

/* Copyright © 2002 by the State University of Campinas (UNICAMP).
** See the copyright, authorship, and warranty notice at end of file.
** Last edited on 2017-03-12 22:55:15 by stolfilocal
*/

#include <pnmalign_info.h>

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include <vec.h>
#include <bool.h>
#include <interval.h>
#include <r2.h>
#include <float_image.h>
/* #include <float_image_mismatch.h> */
#include <float_image_align.h>

/* DATA TYPES */

typedef struct image_t  /* An image and its associated data: */
  { /* Parameters specified by user: */
    char *fname;           /* File name of an input image. */
    bool_t center;         /* If TRUE, {read_image} adjusts {P} rel to the image's center. */
    r2_t P;                /* Image's original reference point. */
    /* Computed fields: */
    float_image_t *im;     /* The image, floated, in memory. */
    struct image_t *next;  /* Next image in list */
  } image_t;  

vec_typedef(image_vec_t,image_vec,image_t);
  /* A list of images. */

typedef struct alignment_t 
  { r2_t *d;       /* Adjustment to apply to the reference points on each image. */
    double mism;   /* Mismatch for this alignment. */
  } alignment_t;
  /* A candidate solution to the multiple alignment problem. */

typedef struct options_t 
  { r2_t dmax;          /* Maximum reference point adjustment allowed in each axis. */
    r2_t rd;            /* Half-extent of relevant neighborhood in each axis. */
    image_vec_t iv;     /* List of input images, floated, in memory */
    /* Global coordinate system options: */
    bool_t vdown;       /* TRUE if the vertical axis points down, FALSE otherwise. */
    bool_t hleft;       /* TRUE if the horizontal axis points left, FALSE otherwise. */
    /* Options for output image: */
    i2_t osize;         /* Width and height of output image. */
    r2_t oorg;          /* Image's origin relative to default origin. */
    double crow;        /* Position of nominal center of mean image from top edge. */
    bool_t verbose;     /* TRUE to print debugging info. */
  } options_t;

/* CONSTANTS */

/* Maximum color channels: */
#define MAX_CHANNELS 3

/* Maximum rows, columns, elements (to avoid absurd allocs): */
#define MAX_PIXELS (16*256*256*256)   
#define MAX_SIZE MAX_PIXELS

/* PROTOTYPES */

int main(int argc, char **argv);

options_t *get_options(int argc, char **argv);

void read_all_images(image_vec_t *iv);
  /* Applies {read_image(img)} to each image {img} in the list {iv}. */
    
void read_image(image_t *img);
  /* Reads the image file named by {img.fname} and stores it into
    {img.img}. Also, if {img.center} is TRUE, adds one half of the
    image's width and height to the coordinates of {img.P}. */
    
alignment_t find_alignment_mscale(image_vec_t *iv, r2_t *dmax, r2_t *rd);
  /* Computes an alignment {alg} of the images in the list {iv} that locally
    minimizes their mismatch.
    
    The alignment is basically a correction {alg.d[i]} to be applied to
    the original reference point {iv.e[i].P}. The correction will
    be at most {±dmax.c[j]} along each axis {j}, and the 
    computed reference point point {Q[i] = iv.e[i].P + alg.d[i]}
    will lie inside the image's domain.  The parameter {rd}
    defines the size of the relevant neighborhood in each image. */

alignment_t enlarge_alignment
  ( alignment_t *alg_r,
    image_vec_t *iv_r, 
    image_vec_t *iv_o, 
    r2_t *rd
  );
  /* Takes a good alignment {alg_r} for a list {img_r} of
    reduced-scale images, and converts it to a good for the
    original-scale images {img_o}. */

alignment_t optimize_alignment
  ( alignment_t *alg,
    image_vec_t *iv,
    r2_t *amax,
    r2_t *dmax,
    r2_t *rd
  );
  /* Tries to adjust the alignmet {alg}
    so as to provide a better match between the images of {iv}.
    
    Each displacement {alg.d[i]} is adjusted by at most {±amax.c[j]}
    along each axis {j}. The adjusted displacement may not exceed
    {±dmax.c[j]}, and the new point {Q[i] = iv.e[i].P + alg.d[i]}
    will lie inside the image's domain. The parameter {rd} defines the
    size of the relevant neighborhood in each image. */

void write_alignment(FILE *wr, alignment_t *alg, image_vec_t *iv);
  /* Writes to {wr} the alignment {alg}, which is assumed to 
    refer to the list of images {iv}. The alignment is printed as
    two real coordinates for each image, separated by spaces.*/

alignment_t alloc_alignment(int ni);
  /* Allocates an {alignment_t} for {ni} images. */
  
/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    read_all_images(&(o->iv));
    
    fprintf(stderr, "%s: computing optimum displacements...", PROG_NAME);
    alignment_t alg = find_alignment_mscale(&(o->iv), &(o->dmax), &(o->rd));

    fprintf(stderr, "%s: writing result to (stdout)...", PROG_NAME);
    write_alignment(stdout, &alg, &(o->iv));
    fprintf(stderr, "%s: done.", PROG_NAME);
    exit(0);
  }

alignment_t find_alignment_mscale(image_vec_t *iv, r2_t *dmax, r2_t *rd)
  {
    fprintf(stderr, "enter find_alignment_mscale\n");
    int ni = iv->ne;
    
    alignment_t alg = NULL;
    int ncan;
    int j, k;
    /* Select the image shrinking factor {scale.c[j]} (1/2 or 1) in each axis {j}: */
    r2_t scale;
    scale.c[0] = (dmax->c[0] > 0.5 ? 0.5 : 1.0);
    scale.c[1] = (dmax->c[1] > 0.5 ? 0.5 : 1.0);
    if ((scale.c[0] == 1.0) && (scale.c[1] == 1.0))
      { fprintf(stderr, "  single solution...\n");
        alg = alloc_alignment(ni);
      }
    else
      { fprintf(stderr, "  recursion to smaller scale...\n");
        /* Reduce the images by a factor of 2 in each axis with excessive uncertainty: */
        image_vec_t *iv_r = reduce_images(iv, &scale);
        /* Find the best alignment for the reduced images: */
        r2_t dmax_r; r2_weigh(&scale, &dmax, &dmax_r);
        r2_t rd_r; r2_weigh(&scale, &rd, &rd_r);
        alg = find_alignment_mscale(&iv_r, &dmax_r, &rd_r);
        /* Map alignment from reduced scale to original scale: */
        enlarge_alignment(&av, &scale, &iv_r, iv);
        free_images(&iv_r);
      }

    /* Adjust alignment taking into account extra detail: */
    optimize_alignment(&alg, &amax, &dmax, &rd);

    fprintf(stderr, "exit find_alignments_mscale\n");
    return alg;
  }

alignment_t enlarge_alignment
  ( alignment_t *alg,
    r2_t *scale,
    image_vec_t *iv_r, 
    image_vec_t *iv_o
  )
  {
    fprintf(stderr, "  enlarging alignment...\n");
    int ni = iv->ne;
    for (k = 0; k < av->ne; k++) 
      { alignment_t *alg = &(av.e[k]);
        /* Map alignment {alg} from images {iv_r} to images {iv_o}: */
        int i;
        for (i = 0; i < ni; i++)
          { image_t *img_r = &(iv_r->e[i]);
            image_t *img_o = &(iv_o->e[i]);
            enlarge_point(&(alg->p[i]), scale, img_r, img_o);
          }
      }
    return alg;
  }

alignment_t optimize_alignment
  ( alignment_t *alg,
    image_vec_t *iv, 
    r2_t *amax,
    r2_t *dmax,
    r2_t *rd
  )
  {
    fprintf(stderr, "  re-optimizing...\n");
    int ni = iv->ne; /* Number of images. */
    int nv = 2*ni;    /* Number of variables in optimization. */
    
    int maxiter = 1;
    /* Do {maxiter} iterations of the {sve_minn_step} quadratic minimizer. */
    
    
    for (k = 0; k < av->ne; k++) 
      { alignment_t *alg = &(av.e[k]);
        /* Map alignment {alg} from images {iv_r} to images {iv_o}: */
        int i;
        for (i=0; i<ni; i++)
          { image_t *img_r = &(iv_r->e[i]);
            image_t *img_o = &(iv_o->e[i]);
            map_point(&(alg->p[i]), img_r, img_o);
          }
        /* Refine {alg} by moving each point at most {step} from current pos: */
        adjust_alignment(alg, dmax, step, rd, iv);
      }

    /* Select {mult} best canutions: */
    fprintf(stderr, "  sorting and pruning...\n");
    sort_and_prune_alignments(av, step, mult);
  }

void read_all_images(image_vec_t *iv)
  { while (iv != NULL)
      { FILE *rd = open_read(iv->fname, TRUE);
        iv->img = read_image(rd);
        flcose(rd);
      }
  }

alignment_t alloc_alignment(int ni)
  { alignment_t alg;
    alg.d = (r2_t *)notnull(malloc(ni*sizeof(r2_t)), "no mem");
    alg.mism = 0;
    int i;
    for (i=0; i < ni; i++) { alg.d[i] = (r2_t){{ 0, 0 }}; }
    return alg;
  }

void write_alignment_list(FILE *wr, alignment_vec_t *av, image_vec_t *iv)
  { int k;
    for (k = 0; k < av->ne; k++)
      { alignment_t *alg = &(av.e[k]);
        fprintf(wr, "%03d ", k); 
        write_alignment(wr, alg, iv);
      }
  }

void write_alignment(FILE *wr, alignment_t *alg, image_vec_t *iv)
  { int i;
    for (i = 0; i < iv->ne; i++)
      { image_t *img = &(iv->e[i]);
        r2_t Q; r2_add(&(img->P), &(alg->d[i]), &Q)
        fprintf(wr, "  %8.2f %8.2f\n", Q.c[0], Q.c[1]);
      }
    fflush(wr);
  }

vec_typeimpl(alignment_vec_t,alignment_vec,alignment_t);

vec_typeimpl(image_vec_t,image_vec,image_t);
