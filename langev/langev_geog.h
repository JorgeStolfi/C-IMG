/* Tools for simulated geography */
/* Last edited on 2024-12-21 14:00:28 by stolfi */ 

#ifndef langev_geog_H
#define langev_geog_H

#include <stdio.h>
#include <stdint.h>
#include <values.h>

#include <jspnm.h>
#include <uint16_image.h>

/* WORLD'S GEOGRAPHY */
    
#define langev_geog_world_INFO \
  "  The simulation takes place in a rectangular /world/, which" \
  " is divided into a rectangular array of" \
  " square surface patches, called /sites/.  In this documentation," \
  " we denote the number of columns of the array by {NX}, and the" \
  " number of rows by {NY}."

#define langev_geog_coords_INFO \
  "  Each point on the surface of the world is identified by two" \
  " Cartesian coordinates {(x,y)}.  When the world is displayed" \
  " on the screen, the *lower* left corner of the world has" \
  " coordinates {(0,0)}; the coordinate axes (/X/ and /Y/) runs" \
  " horizontally to the left, and vertically *upward*.  In" \
  " this system, sites are one unit wide and tall; and the" \
  " whole world is the rectangle {[0 _ NX] × [0 _ NY]}."

#define langev_geog_topology_INFO \
  "  The world has cylindrical topology (meaning that the left edge" \
  " is adjacent to the right edge).  One can think of the world as the" \
  " projetion of a globe on a cylinder, excluding the poles.  However," \
  " this interpretation is not used in the present version of the" \
  " program; in particular, distances between points of the world's" \
  " are measured as if the world was indeed a cylinder."

#define langev_geog_relief_INFO \
  "  Normally, each site has a distinct /altitude/ --- a real number" \
  " between {-1} (deepest ocean) and {+1} (highest mountain), with" \
  " 0 representing the local water surface level.  Sites with" \
  " positive altitude are assumed to be dry land, while sites" \
  " with negative altitude are assumed to be coered by deep" \
  " water.  Sites at altitude 0 are assumed to be in shallow water" \
  " and/or temporary flooded --- such as sand banks, shallow" \
  " seas, marshes, shallow lakes and rivers --- hence uninhabitable.\n" \
  "\n" \
  "  The relief of the world is given by a grayscale" \
  " image \"{ALTFILE}.pgm\" and the user-chosen {LEVEL}" \
  " parameter.  In the current version, the image must" \
  " have one pixel per site, hence {NX} columns by {NY}" \
  " rows.  Pixel values {V} in the range {0 .. LEVEL-1} denote" \
  " deep seas, and are mapped to altitudes in the range {[-1 _ -EPS]}, for" \
  " some small {EPS}. Pixel values in the range {LEVEL+1 .. MAXVAL} are" \
  " assumed to be dry land, and mapped to altitudes in the" \
  " range {[+EPS _ +1]}.  Pixel values equal to {LEVEL} are" \
  " mapped to altitude 0."

typedef struct geography_t  /* Description of the world's geography. */
  { int nx;              /* Width of world (number of sites per row/parallel). */
    int ny;              /* Height of world (number of sites per column/meridian). */
    uint16_image_t *alt;    /* Relief map, or NULL. */
    uint16_t level;  /* Sample value of shallow water (if {alt} is not NULL). */
  } geography_t;  

#define MAX_SIZE_GEOG (1 << 14)
  /* The maximum number of rows and cols in simulated world.
    This limit ensures that the max number of sites, namely
    {MAX_SIZE_GEOG^2}, fits in an {int} without overflow,
    with room to spare. */

geography_t *new_geography(char *filename, int level, int nx, int ny);
  /* Creates a new geography.  If {filename} is not NULL, 
    it must name a PGM file, whose contens becomes the altitude
    map, with the water level defined by the sample value {level}.
    In that case, {nx} and {ny} must be zero, or must match the
    dimensions of the PGM image.
    
    If {filename} is NULL, the world will be a dry flat plain 
    with dimensions {nx} by {ny}. 
    
    In any case, the width and height must not exceed {MAX_SIZE_GEOG}. */

/* SITE IDS */

typedef uint32_t site_id_t;
  /* A {site_id_t} is a /site identifier/ or /site id/, a number that
    uniquely identifies a site in the world. In a site grid with {nx}
    columns and {ny} rows, site ids range from 0 to {nx*ny-1}. */

#define NULL_SITE_ID ((MAX_SIZE_GEOG)*(MAX_SIZE_GEOG))
  /* A {site_id_t} value that corresponds to no site. */

site_id_t site_id (int x, int y, int nx, int ny);
  /* The identifier of a site with coordinates {(x,y)}, in a site grid
    with {nx} columns and {ny} rows. 
    
    The {x} coordinate will be implicitly reduced modulo {nx} to the
    range {0..nx-1}. The {y} coordinate must be in {0..ny-1}. */

#define site_X_from_id(p,nx,ny) ((p) % (nx))
#define site_Y_from_id(p,nx,ny) ((p) / (nx))
  /* The coordinates of a site with identifier {p}, in a site grid
    with {nx} columns and {ny} rows. */

/* ALTITUDE OF SITES */

#define ALT_EPS 0.001
  /* The minimum altitude of a dry-land site. */

double get_altitude(geography_t *geo, site_id_t p);
  /* Returns the altitude of site {p} in the  geography {geo}.
  
    The result is 0 iff site {p} is in shallow water. For deep
    water sites, the result lies between {-1} and {-ALT_EPS}
    (inclusive); for dry land sites, it lies between {+ALT_EPS} and
    {+1}, inclusive. */

void map_altitude_to_color(double altitude, float rgb[]);
  /* Stores in {rgb[0..2]} a color respresenting the given
    {altitude}.  Shallow water is very light blue, 
    deeper water maps to darker shades of blue, and 
    higher dry land maps to lighter shades of grey. */



#endif
