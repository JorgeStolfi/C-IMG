/* Last edited on 2024-12-20 21:01:53 by stolfi */

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

        
        
          }
        else
          { /* This chunk overlaps the previous one. */
            assert(nchunk_old >= 0);
            assert(ichunk_old >= 0);
            assert(ichunk >= nread - nchunk_old);
            /* Shift down the overlapped portion in {si}: */
            uint32_t nover = nread - ichunk;
            uint32_t nskip = nchunk_old - nover;
            fprintf(stderr, "reusing %d samples [%d..%d]\n", nover, nread - nover, nread - 1);
            for (uint32_t k = 0; k < nover; k++)
              { for (uint32_t c = 0; c < nc; c++)
                  { si.sv[c][k] = si.sv[c][nskip + k]; }
              }
            /* Read the new portion: */
            uint32_t nrest = nchunk - nover;
            fprintf(stderr, "reading %d samples [%d..%d]\n", nrest, nread, nread + nrest - 1);
            jsa_read_au_file_samples(stdin, &h, &si, nover, nrest);
            nread += nrest;
          }
