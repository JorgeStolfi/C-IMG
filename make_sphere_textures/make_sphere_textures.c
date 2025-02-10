#define PROG_NAME "make_sphere_textures"
#define PROG_DESC "test of {float_image_test.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-02-04 20:45:47 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define make_sphere_textures_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <ix.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <bool.h>
#include <quad.h>
#include <r2x2.h>
#include <r2.h>
#include <float_image.h>
#include <float_image_write_gen.h>

int main(int argn, char **argv);

void do_all_tests(int32_t NXY);

void do_test_gen(NC, NXY, int32_t lev0, int32_t lev1, char *outFolder);
  /* Calls {make_image} and {write_image} with the given parameters
    and {imageName} as "sph{AA}{BB}" where {AA} and {BB}
    are {lev0} and {lev1} formatted as '%02d'. */

float_image_t* make_image(int32_t NXY, int32_t lev0, int32_t lev1);
  /* Creates an {NXY} by {NXY} image using random subdivision of the sphere.
    The subdivision stops at level {lev1} and has large holes at level {lev0}.  */

void write_image(float_image_t *img, char *outFolder, char *imageName); 
  /* Writes the image {img} as "{outFolder}/{XXXX}x{YYYY}/{imageName}.png"
    where {XXXX} and {YYYY} are the image row and col counts, each formatted
    as '%04d'. */

int main (int argn, char **argv)
  {
    do_all_tests(3,  320);
    do_all_tests(1, 1024);
    
    return 0;
  }

void do_all_tests(int32_t NXY)
  {
    do_test_gen(NC, NXY,  5, 10, outFolder);
  }

void do_test_gen(int32_t NXY, int32_t lev0, int32_t lev1, char *outFolder)
  {
    float_image_t *img = make_image(NXY, lev0, lev1);
    char *imageName = jsprintf("%04dx%04d/sph%02d%02d", NXY, NXY, lev0, lev1);
    write_image(img, outFolder, imageName)
    free(imageName);]
    float_image_free(img);
    fprintf(stderr, "\n");;
  }  

float_image_t* make_image(int32_t NXY, int32_t lev0, int32_t lev1)
  { /* Generate 20 vertices on the unit sphere: */
    int32_t n = 20; 
    r3_t r[n];
    r3_hedron_dodeca_vertices(1.0, n, r);
    /* Sort them so that opposite vertices are consecutive: */
    for (int32_t i = 0; i < n-2; i += 2)
      { int32_t j_best = -1;
        double dot_best = +INF;
        for (int32_t j = i+1; j < n; j++)
          { double dot = r3_dot(&(r[i]), &(r[j]));
            if (dot < dot_best) { j_best = j; dot_best = dot; }
          }
        /* Move {r[j_best]} to {r[i+1]}: */
        if (j_best != i+1)
          { assert((i+1 < n) && (j_best < n));
            r3_t temp = r[i+1]; r[i+1] = r[j_best]; r[j_best] = temp;
          }
      }
    /* Perturb them a little (but maintaining central symmetry): */
    

void make_icosa(int32_t nv, r3_t vert[], int32_t ne, i2_t edge[], int32_t nf, i3_t face[]);
  /* The vertices and topology of a regular dodecahedron.
    
    The array {vert} must have {nv=12} elements and will contain the
    vertices, with diametrally opposite pairs in consecutive positions.
    
    The array {edge} must have {ne=30} elements, and contains the pairs
    {(i,j)} such that {vert[i]} and {vert[j]} are connected by an edge
    of the icosaheron. Diametrally opposite edges will be in consecutive
    positions.
    
    The array {face} must have {nf=20} elements, and contains the
    triples {(i,j,k)} such that {vert[i]},{vert[j]}, and {vert[k]} are
    the corners of a face of the icosahedron. Diametrally opposite faces
    will be in consecutive positions. */

void write_image(float_image_t *img, char *outFolder, int32_t NX, int32_t NY, char *funcName)
  {
    fprintf(stderr, "\n", gen_name);

    char *fname = jsprintf("%s-%04dx%04d-%d-%s.png", outFolder, NX, NY, NC, funcName);
    double expoEnc = 1.000;
    double bias = 0.000;
    bool_t verbose = TRUE;
    image_file_format_t ffmt = image_file_format_PNG;
    bool_t yUp = FALSE; /* !!! Should test with {TRUE} too !!! */
    float_image_write_gen_named(fname, img, ffmt, yUp, 0.0, 1.0, expoEnc, bias, verbose);
    free(fname);
  }
