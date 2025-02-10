#ifndef multifok_make_stack_read_image_H
#define multifok_make_stack_read_image_H
/* Last edited on 2025-01-30 05:04:31 by stolfi */

#include <stdint.h>

#include <bool.h>
#include <float_image.h>
#include <image_file_format.h>

float_image_t *multifok_make_stack_read_image
  ( char *fname, 
    image_file_format_t ffmt, 
    bool_t yUp, 
    bool_t verbose
  );
  /* Reads the reference image file whose name is {fname}, converting it
    from the format {ffmt} to the in-memory
    {float_image_t} format, as described in {float_image_read_gen_INFO}
    and {float_image_read_gen_CONV_INFO}.  */

#endif
