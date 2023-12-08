/* K.Glasov's FFT implementation. */
/* Last edited on 2023-03-16 23:02:38 by stolfi */

#ifndef odt_fft_H
#define odt_fft_H

#define _GNU_SOURCE
#include <stdint.h>

#include <odt_basic.h>

int32_t Rev(int32_t n, int32_t x);
  /* Reverses the low-order {n} bits of {x}. */

res RevSort (int32_t n, cmp * (wherefrom), cmp * (whereto));
  /* Rearranges {wherefrom[0..n-1]} into {whereto[0..n-1]} according to
    the bit-reversal permutation. */

void FFT2d (int32_t way, int32_t m, int32_t n, cmp * ax, cmp * ay);
  /* FFT of a 2D array. */
  
void TR1 (int32_t m, int32_t n, cmp * ax, cmp * ay);
  /* Stores into {ay} the transpose of the matrix {ax}. */

void TR (int32_t n, int32_t m, cmp * ax);
  /* Transposes the matrix {ax}.  Assumes square? */

void FFT (int32_t way, int32_t n, cmp * ax, cmp * ay);
  /* Places in {ay[0..n-1]} the FFT of {ax[0..n-1]}. */

void Lanco (int32_t way, int32_t n, cmp * wherefrom, cmp * whereto);
  /* Stage {way} of the Fourier transform. */

cmp WWW (int32_t way, int32_t n, int32_t k);
  /* Returns {w^(way*k)} where {w} is the primitive {n}th root of uinty. */

int32_t log2i (int32_t x);
  /* Given a positive integer {x}, returns the {log} base 2 of {x}, 
    rounded down: that is, the number of bits needed to store {0..x}. */

#endif
