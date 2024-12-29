/* Utilities */
/* Last edited on 2024-12-21 14:00:04 by stolfi */ 

#ifndef langev_util_H
#define langev_util_H

#include <stdio.h>
#include <stdint.h>
#include <values.h>

#include <bool.h>
#include <vec.h>
#include <uint16_image.h>

#include <langev_lang.h>
#include <langev_gene.h>
#include <langev_name.h>
#include <langev_base.h>
#include <langev_move.h>

#define TRUNCATED_EXP_MIN_Z (-10)

double truncated_exp(double z);
  /* Returns {exp(z)}, slightly tweaked so that it is zero when {z}
    is less than {TRUNCATED_EXP_MIN_Z}. */

#define TRUNCATED_BELL_MAX_Z (4)

double truncated_bell(double z);
  /* Returns the Gaussian-like function {exp(-ln(2)*z^2) == (1/2)^{z^2}}.
    It is 1 when {z==0}, {1/2} when {z==1}, and {1/16} when {z==2}.
    
    The function is slightly tweaked so that it is zero when {|z|} is
    greater than {TRUNCATED_BELL_MAX_Z}. */

#define TRUNCATED_SIGMOID_MAX_Z (4)

double truncated_sigmoid(double z);
  /* Returns the sigmoid function {(1 + erf(z))/2}, that monotonically
    grows from 0 to 1 when {z} goes from {-INF} to {+INF}. At {z==0},
    the function has value {1/2}, slope {2/sqrt(pi)}, and zero
    acceleration.
    
    The function is slightly tweaked so that it is exactly 0 or 1 when
    {|z|} is greater than {TRUNCATED_SIGMOID_MAX_Z}. */

void pick_colors(int M, int32_t x[], int n, float rgb[]);
  /* Picks {n} colors {rgb[3*i+k]}, for {i} in {0..n-1} and {k} in
    {0..2}, so that if {x[i]=x[j]} then color {i} is equal to color
    [j]. Assumes that only the lowest {M} bits of {x[i]} are
    meaningful. Tries to preserve distances while spreading the points
    as best as possible. Does not use gray colors. */

int pick_elem(double w[], int n, int idef, double wtot);
  /* Given a vector {w[0..n-1]} of non-negative weights, returns an
    integer {k} in the range {0..n-1} with probability proportional to
    {w[k]}. 
    
    The procedure returns {idef} if {n} is zero, or if the sum of
    {w[0..n-1]} is zero, or if the selection failed because of
    roundoff errors.
    
    If the argument {wtot} is positive, it must be the sum of
    {w[0..n-1]}.  This option is interesting when the sum is known
    beforehand (e.g. {wtot==1}), or when doing many calls to
    {pick_elem} with the same weight vector. */

typedef void sibs_color_func_t(site_id_t pc, sibs_t *sc, float rgb[]);
  /* Type of any function that stores in {rgb[0..3]} a color value corresponding
    the site {pc}, assuming that it has contents {*sc}. */

uint16_image_t *colorize_frame(frame_t *frm, sibs_color_func_t *color);

void print_hex_alpha(FILE *wr, uint32_t w, int nbits);
  /* Prints the last {nbits} bits of word {w} in hexadecimal 
     using the digits 'A' through 'P' for 0 to 15.
     The digits are printed in order of decreasing significance.
     The word {w} is padded with 0 at the most significant end
     to complete an integer number of digits.  If {nbits} is
     zero, nothing gets printed. */

void print_binary(FILE *wr, uint32_t w, int nbits);
  /* Prints the last {nbits} bits of word {w} in binary.
    The bits are printed in order of decreasing significance.
    If {nbits} is zero, nothing gets printed. */

#endif
