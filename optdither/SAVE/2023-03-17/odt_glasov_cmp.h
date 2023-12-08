#ifndef odt_glasov_cmp_H
#define odt_glasov_cmp_H

/* Complex artithmetic tools by Konstantin Glasov, 1993. */
/* Last edited on 2023-03-17 06:02:08 by stolfi */ 

/* COMPLEX MATH */

#define _GNU_SOURCE
#include <stdint.h>

typedef struct cmp { double re; double im; } cmp;
  /* A complex number in standard form. */

typedef struct plr { double r; double an; } plr ;
  /* A complex number in polar form. */

cmp A(cmp x, cmp y);
  /* Returns the sum {x+y}. */

cmp M(cmp x, cmp y);
  /* Returns the product {x+y}. */

cmp mPolar(plr x);
  /* Converts the complex {x} from polar form to standard form. */

double mAbs(cmp x);
  /* Absolute value (modulus) of {x}. */

plr mExp_p(cmp x);
  /* ???. */

plr mCmp(cmp x);
  /* Converts the complex number {x} from standard form to polar form. */

double mCmpdd(int32_t m,int32_t n,cmp *wfrom,plr *wto);
  /* Converts the array of standard complexes {wfrom[0..m*n-1]} to an
    array {wto[0..m*n-1]} of polar forms, and returns the maximum
    modulus of those numbers. */

double Power(int32_t m, int32_t n, cmp *x);
  /* Converts the array of standard complexes {wfrom[0..m*n-1]} to an
    array {wto[0..m*n-1]} of polar forms, and returns its total
    power. */

#endif
