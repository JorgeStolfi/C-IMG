#define PROG_NAME "imqtopgm"
#define PROG_DESC "convert NASA-JPL-PDS's \".imq\" image files to PGM"
#define PROG_VERS "1.1"

/* Last edited on 2024-12-21 14:00:46 by stolfi */

/* Copyright © 1996 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */

#define PROG_HELP \
  PROG_NAME " < INFILE.imq > OUTFILE.pgm"

#define PROG_INFO \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads from {stdin} an image in the NASA/JPL/PDS \".IMQ\" format, used" \
  " by the Voyager and Viking probes. Writes to {stdout} the image in" \
  " raw PNM grayscale (\".pgm\") format, with 8 bits per pixel.\n" \
  "\n" \
  "  The conversion entails decompressing and integrating the Huffman-encoded" \
  " pixel differences along each line.  Other than that, the pixels are not" \
  " gamma-adjusted or processed in any way.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 1996-08-04 by Jorge Stolfi, IC-UNICAMP,\n" \
  "   based on \"vdcomp.c\" by Mike Martin for NASA-JPL-PDS, 1988-1989.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  1988-07 M.Martin: Wrote \"vdcomp.c\" to decompress Voyager images.\n" \
  "  1989-06 M.Martin: Fixed READVAR to get length on 16-bit unswapped hosts.\n" \
  "  1989-08 M.Martin: Added filenames as arguments, freeing Huffman tree memory.\n" \
  "  1989-12 M.Martin: Handling of Viking images as well as Voyager.\n" \
  "  1996-08 J.Stolfi: Radically simplified the code; output \".pgm\"only." \
  "  2023-02 J.Stolfi: Complete rewrite using JSLIBS libraries." \
  "  2023-02 J.Stolfi: Replaced Kris Becker's \"DECOMP.C\" by \"huff_tree.h\".\n" \
  "\n" \
  "  See the \"libjs\" library in the JSLIBS package on GitHub." \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <argparser.h>
#include <codetree.h>
#include <huff_tree.h>
#include <imq_huff.h>

#include <imq.h>

typedef uint8_t byte_t;

/* Procedure prototypes: */


void write_pgm_header
  ( FILE *wr, 
    uint32_t NX, 
    uint32_t NY,
    bool_t verbose
  );
  /* Writes to {wr} the header of a PGM "raw" (P5) image file
    with {NX} columns and {NY} lines. */

uint32_t small_endian_uint32(byte_t b[]);
  /* Parses {b[0..4]} as a {unit32_t} value in small-endian format. */
    
void skip_image_histogram(FILE *rd, char format, bool_t verbose);
  /* Reads from {rd} the image histogram. In Voyager images it is 2 records,
    in the Viking images it is 1 record. The histogram is currently ignored. */
  
void read_pixel_diff_histogram(FILE *rd, char format, uint32_t nh, uint64_t hist[], bool_t verbose);
  /* Reads the histogram {hist[0..nh-1]} of pixel differences used to build the 
    Huffman tree decoder.
  
    Pixel differences are integers in {-255..+255}, thus {nh} should be
    511. The frequency of the difference {d} will be {hist[d+255]}. */
    
void skip_eng_summary(FILE *rd, bool_t verbose);
  /* Skips the engineering summary record. */
      
void skip_line_header_table(FILE *rd, bool_t verbose);
  /* Skips the line header table (Viking images only). */
      
bool_t read_and_decompress_image_line(FILE *rd, uint32_t line, codetree_t *tree, uint32_t np, byte_t pix[], bool_t verbose);
  /* Reads a compressed image line, decompresses it with the given {tree},
    integrates the pixel differences, and stores the pixels into {pix[0..np-1]}.
    Returns {TRUE} if it hits end-of-file while trying to read the next record.
    
    Fails if the line does not have 1 explicit pixel (byte) plus {np-1}
    Huffman-encoded differences. */

void write_pgm_line(FILE *wr, uint32_t np, byte_t pix[], bool_t verbose);
  /* Writes the image samples {pix[0..np-1]} to {wr}, in the "raw" PGM format. */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    bool_t verbose = TRUE;
    bool_t debug = FALSE;

    FILE *rd = stdin;
    FILE *wr = stdout;

    char format; /* File format: 'Y' for Voyager, 'K' for Viking. */
    uint32_t record_bytes;
    uint32_t NX; /* Image width (number of 8-bit pixels per line). */
    uint32_t NY; /* Image height (number of lines). */
    uint32_t label_checksum;
    imq_read_header(rd, &format, &record_bytes, &label_checksum, &NX, &NY, verbose);

    if (verbose) { fprintf(stderr, "writing output PGM (raw) header...\n"); }
    write_pgm_header(wr, NX, NY, verbose);

    skip_image_histogram(rd, format, verbose);
    
    uint32_t nh = 511;
    uint64_t hist[nh]; /* Histogram of pixel differences as {uint64_t} counts. */
    read_pixel_diff_histogram(rd, format, nh, hist, verbose);

    skip_eng_summary(rd, verbose);

    if (format == 'K') { skip_line_header_table(rd, verbose); }

    /* Initialize the decompression: */
    codetree_t *tree = imq_huff_build_tree(hist);
    if (debug) 
      { fprintf(stderr, "imqtopgm: bit codes:\n");
        imq_huff_print_codes(stderr, tree);
      }

    /* Decompress the image, line by line: */
    uint32_t checksum = 0;
    uint32_t line = 1;
    byte_t pix[record_bytes];
    while(line <= NY)
      { uint32_t np = (uint32_t)NX;
        bool_t verbline = debug && (line == 1);
        if (verbline) { fprintf(stderr, "imqtopgm: decompressing line %d\n", line); }
        bool_t eof = read_and_decompress_image_line(rd, line, tree, np, pix, verbline);
        if (eof) { break; }
        line += 1;
        write_pgm_line(wr, np, pix, verbline);
        if (format == 'K') /* do checksum for viking */
          { for (uint32_t i = 0;  i < np; i++)
              { checksum += (uint32_t)pix[i]; }
          }
      }
    if ((format == 'K') && (label_checksum != checksum))
      { fprintf(stderr, "imqtopgm: !! warning: image checksum expected = %u computed = %u\n", label_checksum, checksum); }
  
    (void)codetree_free(tree);
    fclose(rd);
    fclose(wr);
    return (0);
  }
   
void skip_image_histogram(FILE *rd, char format, bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "imqtopgm: skipping image histogram...\n"); }
    uint32_t nb_max = 4*256; 
    byte_t bbuf[nb_max]; /* Histogram as bytes read from file, small-endian 4-byte nums. */
    uint32_t nb_read = 0;
    uint32_t length1 = 0, length2 = 0;
    length1 = imq_read_var_length_record(rd, nb_max, bbuf);
    if (verbose) { fprintf(stderr, "imqtopgm: read one record (%u bytes)\n", length1); }
    nb_read = nb_read + length1;
    if (format == 'Y') /* read one more time for Voyager image */
      { /* In Voyager images the image histogram is two records of 836 and 188 bytes: */
        demand(length1 == 836, "invalid Voyager record length (image hist 1)"); 
        length2 = imq_read_var_length_record(rd, nb_max-nb_read, bbuf+nb_read);
        if (verbose) { fprintf(stderr, "imqtopgm: read another record (%u bytes)\n", length2); }
        demand(length2 == 188, "invalid Voyager record length (image hist 2)");
        nb_read = nb_read + length2;
        
      }
    else if (format == 'K')
      { /* In Viking images, the image histogram is a single record of 1024 bytes: */
        demand(length1 == 1024, "invalid Viking record length (image hist)");
      }
    else
      { assert(FALSE); }
      
    if (verbose) 
      { fprintf(stderr, "imqtopgm: skipped total %d + %d = %d bytes", length1, length2, nb_read);
        fprintf(stderr, " (%d entries)", nb_read/4); 
        fprintf(stderr, " of image histogram\n"); 
      }
  }
  
void read_pixel_diff_histogram(FILE *rd, char format, uint32_t nh, uint64_t hist[], bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "imqtopgm: reading the pixel difference histogram...\n"); }
    demand(nh == 511, "wrong histogram size");
    uint32_t nb_max = 4*nh;
    byte_t bbuf[nb_max]; /* Histogram as bytes read from file, small-endian 4-byte nums. */
    uint32_t length1 = 0, length2 = 0, length3 = 0;
    uint32_t nb_read = 0;
    length1 = imq_read_var_length_record(rd, nb_max-nb_read, bbuf+nb_read); nb_read += length1;
    if (verbose) { fprintf(stderr, "imqtopgm: read one record (%u bytes)\n", length1); }
    length2 = imq_read_var_length_record(rd, nb_max-nb_read, bbuf+nb_read); nb_read += length2;
    if (verbose) { fprintf(stderr, "imqtopgm: read another record (%u bytes)\n", length2); }
    if (format == 'Y')
      { /* In Voyager images, the difference histogram is three records of 836, 836, and 372 bytes: */
        demand(length1 == 836, "invalid Voyager record length (diff hist 1)"); 
        demand(length2 == 836, "invalid Voyager record length (diff hist 2)"); 
        length3 = imq_read_var_length_record(rd, nb_max-nb_read, bbuf+nb_read); nb_read += length3;
        if (verbose) { fprintf(stderr, "imqtopgm: read another record (%u bytes)\n", length3); }
        demand(length3 == 372, "invalid Voyager record length (diff hist 3)"); 
      }
    else if (format ==  'K')
      { /* In Viking images, the image histogram is two records of 1024 and 1020 bytes: */
        demand(length1 == 1204, "invalid Viking record length (diff hist 1)");
        demand(length2 ==  840, "invalid Viking record length (diff hist 2)");
      }
    else
      { assert (FALSE); }
    if (verbose) 
      { fprintf(stderr, "imqtopgm: read diff histogram"); 
        fprintf(stderr, " as %d + %d + %d = %d bytes", length1, length2, length3, nb_read); 
        fprintf(stderr, " (%d values)\n", nb_read/4);
      }
    demand(nb_read == 4*nh, "diff histogram length is not 511");

    /* Parse {bbuf} into 4-byte integers: */
    for (uint32_t i = 0;  i < nh; i++) 
      { hist[i] = (uint64_t)small_endian_uint32(&(bbuf[4*i]));
        if (verbose) { fprintf(stderr, "  %+4d %lu\n", i-255, hist[i]); }
      }
      
  }

void skip_eng_summary(FILE *rd, bool_t verbose)  
  { uint32_t nb_max = 2048;
    byte_t ibuf[nb_max];
    uint32_t length = imq_read_var_length_record(rd, nb_max, ibuf);
    if (verbose) { fprintf(stderr, "imqtopgm: skipped %d bytes of engineering summary\n", length); }
  }

void skip_line_header_table(FILE *rd, bool_t verbose)
  { uint32_t nb = 2048; /* Max bytes per record. */
    byte_t ibuf[nb];
    uint32_t nr = 1056; /* Number of records to skip. */
    uint32_t tot_bytes = 0;
    for (uint32_t i = 0;  i < nr; i++)
      { uint32_t length = imq_read_var_length_record(rd, nb, ibuf);
        tot_bytes += length;
      }
    if (verbose) { fprintf(stderr, "imqtopgm: skipped %u bytes in %d records of line header table\n", tot_bytes, nr); }
  }

bool_t read_and_decompress_image_line(FILE *rd, uint32_t line, codetree_t *tree, uint32_t np, byte_t pix[], bool_t verbose) 
  { 
  
    uint32_t nb_max = 2048;
    byte_t ibuf[nb_max];
    uint32_t length = imq_read_var_length_record(rd, nb_max, ibuf);
    if (length <= 0) { return TRUE; }
    if (verbose) 
      { fprintf(stderr, "imqtopgm: decompressing line %u", line);
        fprintf(stderr, "  read %u bytes\n", length);
        fprintf(stderr, "  expecting %u pixels...\n", np); 
      }
    
    /* The first byte of {ibuf} is taken to be {pix[0]}. The remaining bytes are interpreted as a
      bit string that is a sequence of differences {diff[0..nd-1]}, Huffman-coded according to 
      the given Huffman code tree. */
  
    if (verbose) { fprintf(stderr, "imqtopgm: decompressing pixels...\n"); }
    imq_huff_decode(line, length, ibuf, tree, np, pix);
    return FALSE;
  }
  
void write_pgm_line(FILE *wr, uint32_t np, byte_t pix[], bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "writing 1 item with %u bytes...\n", np); }
    uint32_t count = (uint32_t)fwrite(pix, np, 1, wr);
    if (verbose) { fprintf(stderr, "wrote %u items...\n", count); }
    demand(count == 1, "error writing output file");
  }


void write_pgm_header
  ( FILE *wr, 
    uint32_t NX, 
    uint32_t NY,
    bool_t verbose
  )
  {
    fprintf(wr, "P5\n");
    fprintf(wr, "%u %u\n", NX, NY);
    fprintf(wr, "255\n");
    fflush(wr);
  }    

uint32_t small_endian_uint32(byte_t b[])
  { uint32_t val = 0;
    for (int32_t i = 3; i >= 0; i--) val = (val << 8) | b[i];
    return val;
  }  
