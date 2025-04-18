#ifndef fvw_texture_H
#define fvw_texture_H

/* fvw_texture.h - tools for colorizing terrains. */
/* Last edited on 2025-01-21 10:13:56 by stolfi */

#include <stdint.h>

#include <frgb.h>
#include <float_image.h>

typedef struct fvw_texture_t
  { float_image_t *timg;
    int32_t chans[3]; 
    bool_t cell;
    float vmin;
    float vmax;
  } fvw_texture_t;
  /* A description of a texture source, consisting of channels
    {chans[0..2]} of image {timg}. See {fvw_texture_get_color}
    for the semantics. */

frgb_t fvw_texture_color_from_value(float v, float cm_vmin, float cm_vmax);
  /* Computes artificial RGB color {clr} for a pixel of value {v},
     whose nominal range is {[cm_vmin _ cm_vmax]}.
     
     The procedure maps {v} to a normalized value {z}, then
     maps {z} to RGB by a color palette.
     
     Specifically, the procedure requires {cm_vmin < cm_vmax}. If {cm_vmin} is
     negative and {cm_vmax} is positive, the mapping from {v} to {z}
     takes {[-cm_vex _ +cm_vex]} to {[-1 _ +1]}, where {cm_vex} is
     {max(|cm_vmin|, |cm_vmax|)} Otherwise, the mapping takes {[cm_vmin
     _ cm_vmax]} to {[0 _ 1]} if {cm_max} is positive, or {[-1 _ 0]} if
     {cm_vmin} is negative. */

frgb_t fvw_texture_get_color
  ( fvw_texture_t *tx,
    int32_t x, int32_t y,
    float v_cur,
    float vmin_cur,
    float vmax_cur
  );
  /* Extracts the samples from row {y} of {tx} (which must not be {NULL}).
    
    Either {tx.chans[0..2]} are all are valid channel indices for
    {tx.timg}, or {chans[0]} is a valid index and {chans[1..2]} are
    {-1}.
    
    If {tx.chans[0..2]} are all valid, samples
    {tx.timg[tx.chans[0..2],x,y} are interpreted as an RGB color.
    In this case, the parameters (v_cur,vmin_cur,vmax_cur}
    and the fields {tx.vmin,tx.vmax} are ignored.
    
    If {tx.chans[0]} is non-negative and {tx.chans[1..2]} are 
    {-1}, sample {tx.timg[tx.chans[0],x,y]] is converted to an RGB value
    using an internal color palette, assuming its range is {[tx.vmin _ tx.vmax]}.
    In this case, the parameters (v_cur,vmin_cur,vmax_cur} are ignored.
    
    If {tx} is null, uses sample {v_cur} instead, converted to RGB as above, 
    assuming its range is {[vmin_cur _ vmax_cur]}.
    
    If the flag {tx.cell} is false, the image {tx.timg} should have the
    same size as the input image {himg} being displayed, and each
    element {tx.timg[*,x,y]} applies to the element {himg[*,x,y]}. If
    {tx.cell} is true, the image {tx.timg} must have one col and one row
    less than {himg}, and each element {tx.timg[*,x,y]} applies to the
    cell of the terrain's domain with corners {(x,y)} and {(x+1,y+1)}. */

#endif
