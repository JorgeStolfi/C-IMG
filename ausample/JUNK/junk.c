/* Last edited on 2006-10-29 01:33:23 by stolfi */

/* ---------------------------------------------------------------------- */

sound_t read_signal_header(FILE *rd);
  /* Reads the header of a sound clip from stream {rd}, which must be 
    in the Sun ".au" audio file format.  Leaves the stream positioned
    so that {read_signal_chunk} will read the first sample. */

void read_signal_chunk(FILE *rd, int skip, int ns, sound_t *s);
  /* Reads a sound clip from stream {rd}, which must be in 
    the Sun ".au" audio file format, with {s.nc} channels 
    and .  The audio file header 
    must have stream must
    be positioned to read a sample.  Skips {skip} samples,
    reads {ns} samples per channel, and stores */
  
/* ---------------------------------------------------------------------- */

sound_t read_signal(FILE *rd)
  { fprintf(stderr, "Reading input file");
    fprintf(stderr, " (Sun \".au\" format)\n");
    return read_au_file(rd);
  }

/* ---------------------------------------------------------------------- */

