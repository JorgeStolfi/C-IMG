#define PROG_NAME "ssitopgm"
#define PROG_DESC "convert SSI (SRI Sun Image) files to PGM files"
#define PROG_VERS "1.0"

/* Copyright © 1994 by the State University of Campinas (UNICAMP).
** See the copyright, authorship, and warranty notice at end of file.
** Last edited on 2006-04-14 10:05:49 by stolfi
*/

#define PROG_HELP \
  PROG_NAME " < SSIFILE > PGM FILE"

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MAXWD 1024
#define MAXHT 1024

/*
** Notice that the images in the JISCT dataset get as wide as 640 elements
** and as tall as 512 elements.
**
** SRI SUN IMAGE (.ssi) files definition adapted from Pascal Fua
**
** bytes 0 -- 7   : file format "srisunim"
** bytes 8        : header type 0
** bytes 9        : data type 0 = int, 1 = float
** byte 10        : bits per pixel
** byte 11        : number of dimensions (1,2,3)
** byte 12 -- 13  : xdim
** byte 14 -- 15  : ydim
** byte 16 -- 17  : zdim
** byte 18 -- 255 : (undefined)
*/

# define HEADER_STRING "srisunim"

# define HEADER_TYPE 0
# define DATA_INT 0
# define DATA_FLOAT 1

/* whole structure should take 256 bytes */
struct  srisunim_header
  { 
    char srisunim[8];
    char header_type;
    char data_type;
    char bits_per_pixel;
    char n_dimensions;
    unsigned char x_dimension_h;
    unsigned char x_dimension_l;
    unsigned char y_dimension_h;
    unsigned char y_dimension_l;
    unsigned char z_dimension_h;
    unsigned char z_dimension_l;
    char padding[238];
  };
    
/* PROTOTYPES */

int read_srisunim_header(
  FILE *fp, 
  int *heightp, 
  int *widthp, 
  char **typep
);

int data_properties(
  char *type, 
  int *bytes_per_pixel_p, 
  int *data_type_p
);

/* END PROTOTYPES */

int main (void)
  { 
    char *type;

    int height;
    int width;
    int nread, nwritten;

    int bytes_per_pixel;
    int data_type;
    
    int row;
    
    char *buffer;

    /* Reads .ssi header: */
    if (read_srisunim_header(stdin, &height, &width, &type))
      { fprintf(stderr, "Bad srisunim file header\n");
	return (-2);
      }
      
    fprintf(stderr, "height = %d   width = %d  type = %s\n", height, width, type);

    if (data_properties(type, &bytes_per_pixel, &data_type) < 0)
      { fprintf(stderr, "unknown srisunim type %s\n", type);
	return (-3);
      }
    
    if ((width > MAXWD) || (height > MAXHT))
      { fprintf(stderr, "image too big, %d rows, %d cols\n", height, width);
	return (-7);
      }
    
    if ((data_type != DATA_INT) || (bytes_per_pixel != 1))
      { fprintf(stderr, "cannnot handle this format yet: %d, %d\n", 
          data_type, bytes_per_pixel
        );
	return (-8);
      }

    /* Writes raw .pgm header: */
    fprintf(stdout, "P5\n%d %d\n255\n", width, height);

    /* Read/write loop: */
    buffer = (char *) malloc (width*bytes_per_pixel);
    if (buffer == NULL) 
      { fprintf(stderr, "malloc failed\n");
	return (-9);
      }
    
    for (row = 0; row < height; row++)
      {
	nread = fread(buffer, bytes_per_pixel, width, stdin);
	if (nread != width)
	  { fprintf(stderr, "eof or read error on row %d after %d pixels\n", row, nread);
	    return (-4);
	  }
          
	nwritten = fwrite(buffer, bytes_per_pixel, width, stdout);
	if (nwritten != width)
	  { fprintf(stderr, "write error on row %d after %d pixels\n", row, nwritten);
	    return (-5);
	  }
      }

    fclose(stdin);
    fclose(stdout);
    return (0);
  }

int read_srisunim_header(
  FILE *fp, 
  int *heightp, 
  int *widthp, 
  char **typep
)
{
  struct srisunim_header header;
  int bytes_per_pixel;

  if (fread((char *) &header, sizeof(header), 1, fp) != 1)
    {
      fprintf(stderr, "read_srisunim_header: can't read header\n");
      return (-1);
    }

  if (memcmp(header.srisunim, HEADER_STRING, sizeof(header.srisunim)) != 0)
    {
      header.header_type = 0;
      fprintf(stderr, "read_srisunim_header: wrong type %s\n", header.srisunim);
      return (-2);
    }
  if (header.header_type != 0)
    {
      fprintf(stderr, "read_srisunim_header: wrong header_type %d\n", 
        header.header_type
      );
      return (-3);
    }
    
  if (header.n_dimensions != 2)
    {
      fprintf(stderr, "read_srisunim_header: only handle 2 dimensions (%d)\n", 
        header.n_dimensions
      );
      return (-4);
    }

  *widthp = header.x_dimension_h * 256 + header.x_dimension_l;
  *heightp = header.y_dimension_h * 256 + header.y_dimension_l;

  bytes_per_pixel = header.bits_per_pixel / 8;
  if ((bytes_per_pixel * 8) != header.bits_per_pixel)
    {
      fprintf(stderr, "read_srisunim_header: only handle multiples of 8 bits %d\n",
	header.bits_per_pixel
      );
      return (-5);
    }
      
  switch (header.data_type)
    {
    case DATA_INT:
      switch (bytes_per_pixel)
	{
	case sizeof(char):
	  *typep = "char";
	  break;
	case sizeof(short):
	  *typep = "short";
	  break;
	case sizeof(long):
	  *typep = "long";
	  break;
	default:
          /* handle the case of distinct ints */
	  if (bytes_per_pixel == sizeof(int))
	    *typep = "int";
	  else
	    {
	      fprintf(stderr, "read_srisunim_header: don't handle %d byte integers\n", 
                bytes_per_pixel
              );
	      return (-6);
	    }
	  break;
	}
      break;
    case DATA_FLOAT:
      switch (bytes_per_pixel)
	{
	case sizeof(float):
	  *typep = "float";
	  break;
	case sizeof(double):
	  *typep = "double";
	  break;
	default:
          fprintf(stderr, "read_srisunim_header: don't handle %d byte floats\n", 
            bytes_per_pixel
          );
          return (-7);
	  break;
	}
      break;
    default:
      fprintf(stderr, "read_srisunim_header: unhandled data type: %d\n",
	header.data_type
      );
      return (-8);
      break;
  }
  return (0);
}

int data_properties(
  char *type, 
  int *bytes_per_pixel_p, 
  int *data_type_p
)
{
  int bytes_per_pixel;
  int data_type;
  int failed = 0;

  if (strcmp(type, "char") == 0)
    {
      bytes_per_pixel = sizeof(char);
      data_type = DATA_INT;
    }
  else if (strcmp(type, "short") == 0)
    {
      bytes_per_pixel = sizeof(short);
      data_type = DATA_INT;
    }
  else if (strcmp(type, "long") == 0)
    {
      bytes_per_pixel = sizeof(long);
      data_type = DATA_INT;
    }
  else if (strcmp(type, "int") == 0)
    {
      bytes_per_pixel = sizeof(int);
      data_type = DATA_INT;
    }
  else if (strcmp(type, "float") == 0)
    {
      bytes_per_pixel = sizeof(float);
      data_type = DATA_FLOAT;
    }
  else if (strcmp(type, "double") == 0)
    {
      bytes_per_pixel = sizeof(double);
      data_type = DATA_FLOAT;
    }
  else
    {
      failed = 1;
      bytes_per_pixel = 0;
      data_type = 0;
    }

  *bytes_per_pixel_p = bytes_per_pixel;
  *data_type_p = data_type;

   return (failed ? -1 : 0);
}


/* COPYRIGHT, AUTHORSHIP, AND WARRANTY NOTICE:
** 
**   Copyright © 1994 by the State University of Campinas (UNICAMP).
**
** Original source unknown; header files maybe from Kurt Konolige's SVM.
** Modified mar/1994 at IC-UNICAMP by Marellus Rosa Macêdo and Jorge Stolfi.
**
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this copyright, authorship, and warranty notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.
** This software is provided "as is", WITHOUT ANY EXPLICIT OR IMPLICIT
** WARRANTIES, not even the implied warranties of merchantibility and
** fitness for a particular purpose. END OF NOTICE.
*/
