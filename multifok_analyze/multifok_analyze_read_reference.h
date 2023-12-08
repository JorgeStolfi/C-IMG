#ifndef multifok_analyze_read_reference_H
#define multifok_analyze_read_reference_H
/* Last edited on 2017-10-25 19:04:11 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <float_image.h>
#include <image_file_format.h>

float_image_t *multifok_analyze_read_reference
  ( char *fname, 
    image_file_format_t ffmt, 
    bool_t verbose
  );
  /* Reads the reference image file whose name is {fname}, converting it
    from the format {ffmt} to the in-memory
    {float_image_t} format, as described in {float_image_read_gen_INFO}
    and {float_image_read_gen_CONV_INFO}.  */

#endif
