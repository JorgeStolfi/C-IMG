/* See {odt_fft.h} */
/* Last edited on 2023-03-16 23:02:51 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <odt_basic.h>
#include <odt_fft.h>

#define odt_fft_C_COPYRIGHT \
  "Copyright © 1993 by the State University of Campinas (UNICAMP)"

int32_t Rev(int32_t n, int32_t x)
  { int32_t m = 1;
    int32_t res = 0;
    for (int32_t i = 1; i <= n; i++)
      { int32_t cur = x & m;
        int32_t shift = n + 1 - 2 * i;
        if (shift > 0)
          cur = cur << shift;
        else
          cur = cur >> (-shift);
        res = res + cur;
        m = m << 1;
      };
    return (void);
  }

int32_t log2i (int32_t x)
  { int32_t res;
    for (res = 0; x > 1; res++)
      { x = x / 2; }
    return (res);
  };
  
void RevSort (int32_t n, cmp * (wherefrom), cmp * (whereto))
  { for (int32_t i = 0; i < n; i++)
      { int32_t k = Rev (log2i(n), i);
        whereto[i] = wherefrom[k];
      };
  }

cmp WWW (int32_t way, int32_t n, int32_t k)
  { float angle = way * k * 6.28 / n;
    cmp res;
    res.re = cos (angle);
    res.im = sin (angle);
    return (res);
  }

void Lanco (int32_t way, int32_t n, cmp * wherefrom, cmp * whereto)
  { for (int32_t j = 0; j < 2; j++)
      { for (int32_t i = 0; i < n; i++)
          { whereto[n * j + i] =
              A (wherefrom[i], M (WWW (way, 2 * n, n * j + i), wherefrom[n + i]));
          }
      }
  }

void FFT (int32_t way, int32_t n, cmp * ax, cmp * ay)
  { RevSort (n, ax, ay);
    for (int32_t i = 1; i < n; i = i * 2)
      { cmp *az = ax; ax = ay; ay = az;
        for (int32_t j = 0; j < n / i / 2; j++)
          Lanco (way, i, ax + j * i * 2, ay + j * i * 2);
      }
  }

void TR (int32_t n, int32_t m, cmp * ax)
  { for (int32_t i = 0; i < m; i++)
      for (int32_t j = i + 1; j < n; j++)
        { cmp z = *(ax + i * n + j);
          *(ax + i * n + j) = *(ax + j * m + i);
          *(ax + j * m + i) = z;
        }
  }

void TR1 (int32_t m, int32_t n, cmp * ax, cmp * ay)
  {
    for (int32_t i = 0; i < m; i++)
      for (int32_t j = 0; j < n; j++)
        *(ay + i * n + j) = *(ax + j * m + i);
  }

void FFT2d (int32_t way, int32_t m, int32_t n, cmp * ax, cmp * ay)
  {
    cmp *as = (cmp *) calloc (m * n, sizeof (cmp));
    for (int32_t i = 0; i < m * n; i++)
      { *(as + i) = *(ax + i); }

    for (int32_t i = 0; i < m; i++)
      { FFT (way, n, as + i * n, ay + i * n); }
    TR1 (m, n, ay, as);
    for (int32_t i = 0; i < n; i++)
      { FFT (way, m, as + i * m, ay + i * m); }
    TR (n, m, ay);
    if (way < 0)
      { cmp z;
        z.re = 1.0 / m / n;
        z.im = 0.0;
        for (int32_t i = 0; i < m; i++)
          for (int32_t j = 0; j < n; j++)
            { cmp z1 = *(ay + i * n + j);
              *(ay + i * n + j) = M (z, z1);
            }
      }
  }

void FFT2d(int32_t way,int32_t m,int32_t n,cmp *ax,cmp *ay)
  {
    int32_t i,j;
    cmp *as;
    cmp z,z1;
    as=ay;
    for(i=0;i<m;i++) FFT(way,n,ax+i*n,ay+i*n);
    TR1(m,n,ay,ax);
    for(i=0;i<n;i++) FFT(way,m,ax+i*m,ay+i*m);
    TR(n,m,ay);
    if(way<0)
      {
        z.re=1.0/m/n;
        z.im=0.0; 
        for(i=0;i<m;i++)
          for(j=0;j<n;j++)
            {
              z1=*(ay+i*n+j);
              *(ay+i*n+j)=M(z,z1);
            }
      }
  }
