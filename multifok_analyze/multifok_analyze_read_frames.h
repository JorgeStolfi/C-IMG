#ifndef multifok_analyze_read_frames_H
#define multifok_analyze_read_frames_H
/* Last edited on 2018-01-03 14:20:54 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <vec.h>
#include <float_image.h>
#include <image_file_format.h>

float_image_t **multifok_analyze_read_frames
  ( char *framePattern, 
    int32_t NF, 
    int32_t frameID[], 
    image_file_format_t ffmt,
    int32_t NC,
    int32_t NX,
    int32_t NY,
    bool_t verbose
  );
  /* Reads {NF} image files, converting them from the format {ffmt} to the in-memory
    {float_image_t} format, as described in {float_image_read_gen_INFO}
    and {float_image_read_gen_CONV_INFO}. 
    
    Assumes that the frames are For each frame index {f} in {0..NF-1}, the name of file with index {f} is obtained
    by formatting {frameID[f]} with the formar specification {framePattern}.
    
    Fails if any frame does not have {NC} channels, {NX} columns,
    and {NY} rows. */

#endif
