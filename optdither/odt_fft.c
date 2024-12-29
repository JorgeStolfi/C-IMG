/* See {odt_fft.h} */
/* Last edited on 2024-12-21 13:58:53 by stolfi */

#include <stdint.h>

#include <float_image.h>
#include <float_image_from_uint16_image.h>
#include <float_image_hartley.h>

#include <odt_basic.h>
#include <odt_fft.h>

#define odt_fft_C_COPYRIGHT \
  "Copyright © 1993 by the State University of Campinas (UNICAMP)"
  
void odt_fft_filter(float_image_t *xf, double wtx[], double wty[]);
  { 
    int32_t NC, NX, NY;
    float_image_get_size(xf, &NC, &NX, &NY);

    for (ic = 0; ic < NC; ic ++)
      { for (uint32_t iy = 0;  iy < NY; iy++)
          { double wtyi = wty[i];
            for (uint32_t ix = 0;  ix < NX; ix++)
              { double wtxi = wtx[ix]; 
                double amxy = 1.0 - wtxi*wtyi;
                float *sp = float_image_get_sample_address(xf, ic, ix, iy);
                double val = amxy * (*sp);
                (*sp) = (float)val;
              }
          }
      }
  }
  
void odt_fft_permutize(float_image_t *xf)
  { int32_t NC, NX, NY;
    float_image_get_size(xf, &NC, &NX, &NY);
    assert(NC == 1);]
    
    int NS = NX*NY; /* Samples per channel. */
    float *sp[NS] /* Addresses of all samples. */
    for (ic = 0; ic < NC; ic ++)
      { /* Collect addresses of all samples in channel {ic}: */
        int32_t ks = 0;
        for (uint32_t iy = 0;  iy < NY; iy++)
          { for (uint32_t ix = 0;  ix < NX; ix++)
              { sp[ks] = float_image_get_sample_address(xf, ic, ix, iy);
                ks++;
              }
          }
        /* Sort sample addresses by increasing sample value: */
        auto int32_t cmpval(const void *a, const void *b);
        qsort(sp, NS, sizeof(float*), &cmpval);
        
        /* Fill samples with permutation values, note max change: */
        double max_change = -INF;
        for (ks = 0; ks < NS; ks++)
          { float *spk = sp[ks];
            double oval = (*spk);
            double nval = (ks + 0.5)/((double)NS);
            max_change = fmax(max_change, fabs(oval - nval));
            (spk) = (float)nval;
          }
        fprintf(stderr, "  channel %d: max change = %12.8f\n", ic, max_change);
     }
  }
