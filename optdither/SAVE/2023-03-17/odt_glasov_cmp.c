/* See odt_glasov_cmp.h */
/* Last edited on 2023-03-17 06:02:25 by stolfi */

#define odt_glasov_cmp_C_COPYRIGHT \
  "Copyright © 1993 by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <bool.h>
#include <jsfile.h>
#include <affirm.h>

#include <odt_glasov_cmp.h>

cmp A(cmp x, cmp y)
  {
    cmp res;
    res.re = x.re + y.re;
    res.im = x.im + y.im;
    return res;
  }
  
cmp M(cmp x, cmp y)
  {
    cmp res;
    res.re = x.re*y.re - x.im*y.im;
    res.im = x.im*y.re + x.re*y.im;
    return res;
  }
  
cmp mPolar(plr x)
  {
    cmp res;
    res.re = cos(x.an)*x.r;
    res.im = sin(x.an)*x.r;
    return res;
  }
  
double mAbs(cmp x)
  {
    return hypot(x.re, x.im);
  }

plr mExp_p(cmp x)
  {
    plr res;
    res.r = exp(x.re);
    res.an = x.im;
    return res;
  }

plr mCmp(cmp x)
  {
    plr res;
    res.r = mAbs(x);
    res.an = (double)atan2(x.im, x.re);
    return res;
  }

float mCmpdd(int32_t m, int32_t n, cmp *wfrom, plr *wto)
  {
    int32_t i,j;
    float maxampl;
    plr z;

    maxampl = 0.0;

    for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
        { z = mCmp(*(wfrom+i*n+j));
          if (z.r > maxampl) maxampl = z.r;
          *(wto+i*n+j) = z;
        }

    return(maxampl);
  }

float Power(int32_t m, int32_t n, cmp *x)
  {
    int32_t i,j;
    float power;
    plr z;

    power=0.0;

    for (i = 0; i < m; i++)
      for (j = 0; j < n; j++)
        { z = mCmp(*(x+i*n+j));
          power = power+z.r; /* !!! should add SQUARES !!! */
        };

    return(power);
  }

