#ifndef odt_tools_H
#define odt_tools_H

/* Tools for {optdither}. */
/* Last edited on 2024-12-21 12:01:01 by stolfi */ 

/* COMPLEX MATH */

#include <stdint.h>

#include <bool.h>
#include <float_image.h>

float_image_t *odt_read_dither_matrix(char *fname);
  /* Reads from file {fname} a grayscale image, assumed to be in 
    the PGM format, and converts it to a float image.  
    
    The image's {maxval} must be equal to the number of pixels
    minus one, and every sample value in {0..maxval}
    must occur exacly once.
    
    The dither matrix is assumed to be in linear lum scale.
    Each sample value {s} is converted to the float value 
    {(s + 0.5)/(maxval+1)}. */


void odt_round_image(int32_t m, int32_t n, cmp *xc, int32_t *xi);
  /* Converts the complex matrix {xc[0..m*n-1]} to an integer
    matrix {xi[0..m*n-1]} by rounding the real part of each
    element. */

void PrT(int32_t m, int32_t n, int32_t x[]);
  /* Prints the integer array {x[0..m*n-1]}. */

void SaveCT(char *name, int32_t m, int32_t n, cmp *x, char mode, int32_t length);
  /* Saves a complex array as a table in CT.dat file. 
    If {mode} is 'p', saves in polar form; if 'c', saves in complex form. */

void Dith(int32_t m, int32_t n, float *x, int32_t *y, int32_t mm, int32_t nn, int32_t *matr);
  /* Dithers array {x} using matrix {matr}, output to array {y}. */

void rebuild(int32_t m, int32_t n, cmp * x, cmp * y);
  /* Restores usual frequency order */

void LiveAgua (int32_t m, int32_t n, cmp * x, int32_t *y);
  /* Restores a dot matrix according to groth of values */

void MskFlt(int32_t m, int32_t n, cmp * ax, cmp * ay, cmp * az);
  /* Filts putting 0 where 0 are in the pattern */

void Justifier (int32_t m, int32_t n, cmp * x, cmp * y, cmp * res);
  /* Multiplies every member after flt to restore total energy */

#endif
