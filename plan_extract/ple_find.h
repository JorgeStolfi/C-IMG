/* ple_find.h -- planning the extraction of standard-size sub-images */
/* Last edited on 2007-12-27 02:28:59 by stolfi */

#ifndef ple_find_H
#define ple_find_H

#define ple_find_H_COPYRIGHT \
  "Copyright © 2007 by the State University of Campinas (UNICAMP)"
  
/* This interface provides functions to plan the extraction of 
  one or more fixed-size sub-images from a given image, by an 
  optimal combination of cropping, padding, and scaling. */

#include <stdio.h>

#include <bool.h>
#include <i2.h>
#include <irange.h>

/* INTERNAL TYPES */
 
typedef enum
  { edge_L = 0, /* Left (min X) edge. */
    edge_B = 1, /* Bottom (min Y) edge. */
    edge_R = 2, /* Right (max X) edge. */
    edge_T = 3  /* Top (max Y) edge. */
  } edge_t;
  /* Index of an edge of a rectangle. Note that bit 0 is the 
    perpendicular axis, and bit 1 specifes min or max. */

typedef struct plan_t 
  { irange_t CBox[2]; /* Domain of {C} relative to {I}. */
    irange_t PBox[2]; /* Domain of {P} relative to {C}. */
    i2_t SSize;       /* Dimensions of scaled-down image {S}. */
    i2_t ESize;       /* Dimensions of extracted images {E[r,s]}. */
    i2_t EPos;        /* Shift of {E[0,0]} relative to {S}. */
    i2_t EStep;       /* Displacements between successive images {E[r,s]}. */
    irange_t tr[2];   /* Ranges of tile indices {r,s}. */
    double score;     /* Score of plan (negative if unacceptable). */
  } plan_t;
  /* An image-extraction plan. The displacement {EPos} is the relative
    position of the lower left corners. All coordinates assume that
    the Y axis points UP. Tile indices along axis {ax} lie in the
    range {tr[ax].end[0]..tr[ax].end[1]}. */

/* PROTOTYPES */
 
/* Maximum rows, columns, elements (to avoid absurd allocs): */
#define ple_MAX_PIXELS (16*256*256*256)   
#define ple_MAX_SIZE ple_MAX_PIXELS

#define INF INFINITY

void ple_eval_plans_for_esize
  ( irange_t TBox[], 
    bool_t paddable[], 
    i2_t *ESize, 
    i2_t *EStep, 
    double maxScale, 
    double maxStretch, 
    double maxPad, 
    double maxCut,
    plan_t *pl_best, 
    int debug
  );
  /* Considers the plans {pl} that yield extracted images
    {E[r,s]} of size {ESize}. 
    
    More precisely, enumerates all plans that entail cropping the
    sub-image {T} described by {TBox} to some image {C}, padding it to
    some image {P}, scaling it to some image {S}, and then extracting
    from {S} one or more images {E[r,s]}. Whenever it finds a plan
    {pl} with better score than {pl_best}, updates {pl_best} with it.
    The client must initialize {pl_best} suitable before the call.
    
    For the image {C}, considers all non-empty sub-images of {T}. For
    the image {P}, considers padding {C} along any edge {e} where {C}
    touches the boundary of {T} and {paddable[e]} is TRUE.

    For the image {S}, considers all interesting dimensions {SX × SY}
    such that the scale factor from {C} to {S} is at most {maxScale},
    and the longest side of {S} is stretched by at a factor not
    exceeding {maxStretch} relative to the size that would preserve
    the aspect of {C}.

    The amont by which an image {E[r,s]} extends outside the
    scaled-down rectangle {T} along either axis should not exceed
    {maxPad} times the size of {E} along that axis. The amount by
    which the union of all {E[r,s]} falls short of {S} along either
    axis must not exceed {maxCut} times the extent of {S} along that
    axis. */

i2_t ple_box_size(irange_t B[]);
  /* Returns the width and height of the box {B}, namely 
    {B[ax].end[1] - B[ax].end[0]} along each axis {ax}. */

bool_t ple_empty_box(irange_t B[]);
  /* TRUE iff the rectangle {CBox} is empty (has zero width or height). */

void ple_debug_box(FILE *wr, char *pre, irange_t B[], char *suf);
void ple_debug_size(FILE *wr, char *pre, i2_t *sz, char *suf);
void ple_debug_pos(FILE *wr, char *pre, i2_t *pos, char *suf);
void ple_debug_ixrange(FILE *wr, char *pre, irange_t *v, char *suf);
  /* These procedures print {pre} and the given object to {wr}, 
    followed by {suf}. */

void ple_debug_plan(FILE *wr, char *pre, irange_t TBox[], plan_t *pl, char *suf);
  /* Writes the {T} domain {TBox} and the plan {pl} to {wr}, for debugging.
    The printout is preceded by {pre} and followed by {suf}. */

#endif
