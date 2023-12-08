/* See langev_geog.h. */
/* Last edited on 2017-06-20 20:30:25 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

#include <affirm.h>
#include <bool.h>
#include <jspnm.h>
#include <jsfile.h>
#include <uint16_image.h>

#include <langev_geog.h>

geography_t *new_geography(char *filename, int level, int nx, int ny)
  {
    geography_t *geo = (geography_t*)notnull(malloc(sizeof(geography_t)), "out of memory");
    if (filename != NULL)
      { /* Read relief image, check {nx,ny}: */
        FILE *rd = open_read(filename, TRUE);
        uint16_image_t *img = uint16_image_read_pnm_file(rd);
        fclose(rd);
        demand(img->chns == 1, "relief must be a PGM image");
        geo->alt = img;
        /* Save {level} and check against {maxval}: */
        demand(level < img->maxval, "world has no dry land");
        geo->level = level;
        /* Set {nx,ny} from image and check against given size: */
        demand((nx <= 0) || (img->cols == nx), "wrong number of columns in map");
        demand((nx <= 0) || (img->rows == ny), "wrong number of rows in map");
        geo->nx = img->cols;
        geo->ny = img->rows;
      }
    else
      { geo->alt = NULL; 
        geo->level = 0;
         /* Set {nx,ny} from given size: */
        geo->nx = nx;
        geo->ny = ny;
      }
    /* Check validity of final size: */
    demand((geo->nx > 0) && (geo->nx <= MAX_SIZE_GEOG), "bad world width");
    demand((geo->nx > 0) && (geo->ny <= MAX_SIZE_GEOG), "bad world height");
    return geo;
  }

site_id_t site_id (int x, int y, int nx, int ny)
  { /* Reduce the {x} coordinate to {0..nx-1}: */
    while(x < 0) { x += nx; }
    while(x >= nx) { x -= nx; }
    /* Check the {y} coordinate: */
    demand((y >= 0) && (y < ny), "Y coordinate is out of bounds");
    /* Compute the sid: */
    return x + nx*y;
  }

double get_altitude(geography_t *geo, site_id_t p)
  { if (geo->alt == NULL)
      { /* No altitude field - flat plains just above sea level: */
        return 0.001;
      }
    else
      { /* Get altitude from altitude field, map to [-1 _ +1]: */
        uint16_image_t *alt = geo->alt;
        uint16_t level = geo->level;
        assert(alt->chns == 1);
        int nx = alt->cols;
        int ny = alt->rows;
        demand(p < nx*ny, "invalid site id");
        int x = site_X_from_id(p, nx, ny);
        int y = site_Y_from_id(p, nx, ny);
        uint16_t v = alt->smp[ny - 1 - y][x];
        assert(v <= alt->maxval);
        if (v < level)
          { /* Below sea level: */
            return ((double)(v-level))/((double)level);
          }
        else if (v > level)
          { /* Above sea level: */
            return ((double)(v-level))/((double)(alt->maxval-level));
          }
        else
          { /* Shallow water: */
            return 0.0;
          }
      }
  }

void map_altitude_to_color(double altitude, float rgb[])
  {
    if (altitude == 0)
      { /* Shallow water - paint light blue: */
        rgb[0] = 0.700; 
        rgb[1] = 0.800; 
        rgb[2] = 1.000;
      }
    else if (altitude < 0)
      { /* Deep water - paint with shades of blue (deeper = darker): */
        double s = 1 + altitude; /* Shalowness. */
        double t = -altitude;    /* Depth. */
        rgb[0] = 0.490*s + 0.070*t; 
        rgb[1] = 0.560*s + 0.080*t; 
        rgb[2] = 1.000;
      }
    else
      { /* Dry land - paint with grey shades (higher = lighter): */
        double s = 1 - altitude; /* Shalowness. */
        double t = altitude;     /* Altitude. */
        double gv = 0.700*s + 0.900*t;
        rgb[0] = gv; rgb[1] = gv; rgb[2] = gv;
      }
  }
