// Basic types for stereo computations
// Last edited on 2015-06-12 19:40:47 by stolfilocal

#ifndef stereo_basics_H
#define stereo_basics_H

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <r3.h>
#include <r3x3.h>

void orthonormalize(r3x3_t *R_p);
  /* Turns {*R_p} into an orthonormal matrix, by applying Gram-Schmidt to 
    successive rows. */

void print_r3_t(FILE *f, char *pref, r3_t *v_p, char *suff);
void print_r3x3_t(FILE *f, char *pref, r3x3_t *R_p, char *suff);

#endif
 
