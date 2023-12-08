#define PROG_NAME "ausample"
#define PROG_DESC "extract sample clips from audio files"
#define PROG_VERS "1.0"

/* Copyright © 2006 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2015-10-18 02:54:10 by stolfilocal */

#define PROG_HELP \
  PROG_NAME "\\\n" \
  "    [ -unit [ sample | second ] ] \\\n" \
  "    -take {CLIP_LENGTH} \\\n" \
  "    -every {CLIP_STRIDE} \\\n" \
  "    [ -splice {SPLICE_LENGTH} ] \\\n" \
  "    < {INPUT_FILE} \\\n" \
  "    > {OUTPUT_FILE} \\\n" \
  "    " argparser_help_info_HELP

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads an audio file, extracts a series of" \
  " segments from it, and outputs their concatenation.\n" \
  "\n" \
  "  The input signal is read from standard input and is converted from" \
  " its external format to a sequence of double-precision floating-point" \
  " samples.  Each segment has length {CLIP_LENGTH}, and successive segments" \
  " start {CLIP_STRIDE} apart.  The series is centered over the length of" \
  " the input signal.\n" \
  "\n" \
  "  If the \"-splice\" option is used, the segments are concatenated" \
  " with a soft-edged cutout mask whose sigmoid part extends" \
  " {SPLICE_LENGTH/2} on either side of the nominal cut.  In any" \
  " case, the extracted segments are concatenated by aligning" \
  " their *nominal* endpoints; so that consecutive segments" \
  " actually overlap by {SPLICE_LENGTH}.  The masks is used" \
  " only between segments; the first segment always starts" \
  " abruptly and the last segment ends abruptly.\n" \
  "\n" \
  "  All times can be specified in seconds or in sampling" \
  " steps, depending on the \"-unit\" option.  In either case," \
  " the specified values are rounded to the nearest multiple of" \
  " the sampling step.  The nominal position of any cut" \
  " is halfway between two consecutive samples.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -unit [ sample | second ] \n" \
  "    Specifies whether the unit for all time parameters is" \
  " one second or one sampling step.  The" \
  " default is \"sample\" \n" \
  "\n" \
  "  -take {CLIP_LENGTH} \n" \
  "    Specifies the duration of the segments to extract.  This" \
  " parameter is required.\n" \
  "\n" \
  "  -every {CLIP_STRIDE}\n" \
  "    Specifies the delay between the starting times of two" \
  " consecutive segments.  This parameter is required.\n" \
  "\n" \
  "  -splice {SPLICE_LENGTH} \n" \
  "    Specifies the length of the soft transition in the cutout" \
  " mask, which is also the length of the soft splices in the output" \
  " file.  The default is 0 (meaning sharp cuts and no overlap).\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  spegram(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created oct/2006 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP)."

#define stringify(x) strngf(x)
#define strngf(x) #x

/* Must define _GNU_SOURCE to get the defintion of {asprinf}. */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <values.h>
#include <assert.h>
#include <math.h>

#include <fftw3.h>

#include <bool.h>
#include <vec.h>
#include <nget.h>
#include <fget.h>
#include <argparser.h>

#include <jsaudio.h>
#include <jsaudio_au.h>
#include <jsrandom.h>

/* DATA TYPES */

typedef enum
  { TUN_SECOND,  /* Unit is second. */
    TUN_SAMPLE   /* UNit is the sampling step. */
  } time_unit_t;
  /* Specifies the unit of time. */

/* Command line options: */
typedef struct options_t
  { time_unit_t unit;     /* The time unit (samples, seconds, etc.). */
    double take;          /* Length of each segment. */
    double every;         /* Stride between segments. */
    double splice;        /* Length of splice. */
  } options_t;

/* PROTOTYPES */

int main(int argc, char **argv);
  /* Main program. */

void write_signal(FILE *wr, sound_t *s);
  /* Writes a sound clip to stream {wr}, in the Sun ".au" audio file format. */

double samples_per_unit(time_unit_t unit, double fsmp);
  /* Computes the number of sampling steps coresponding to the {unit},
    assuming that the sampling frequency is {fsmp} (in Hz). */

void cut_and_paste_segment
  ( sound_t *si, 
    int ipos, 
    int ntake, 
    double rlo,
    double rhi,
    sound_t *so, 
    int opos
  );
  /* Extracts a segment from {si} and pastes it onto {so}.
    Nominally, the input segment has {ntake} samples, starts at sample
    {si[ipos]}, and is pasted starting at sample {so[opos]}. 
    
    However, if {rlo} is positive, the cut is smoothed by a Hann step
    of half-width {rlo}. In that case, the operation may affect
    several samples before {si[ipos]} and {so[opos]}. The same thing
    holds for the other end, if {rhi} is positive. Assumes that all
    samples affected by the cut/paste exist. */
 
double smooth_step(double x, double r);
  /* Computes a C1-smooth step function that is 0 for {x <= -r},
    0.5 for {x == 0}, and 1 for {x >= +r}. */

options_t* get_options(int argc, char **argv);
  /* Parses the command-line arguments. */

/* IMPLEMENTATIONS */

int main(int argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    /* Read the input file header: */
    au_file_header_t h = jsa_read_au_file_header(stdin);
    
    /* Compute number of channels {nc} and samples per channel {ns}: */
    int nc = h.channels;
    int bps = jsa_au_file_bytes_per_sample(h.encoding);
    int ns = h.data_size / nc / bps;
    assert(h.data_size == ns * nc * bps); 
    fprintf(stderr, "input channels = %d\n", nc);
    fprintf(stderr, "input samples %d\n", ns);
    
    /* Get sampling frequency: */
    double fsmp = (double)h.sample_rate;
    fprintf(stderr, "sampling frequency = %.4f  Hz\n", fsmp);
    
    /* Get number of samples per unit: */
    double spu = samples_per_unit(o->unit, fsmp);

    /* Compute nominal number {ntake} of samples per segment: */
    int ntake = (int)floor(o->take * spu + 0.5); 
    assert(ntake >= 0);
    double ttake = ((double)ntake)/fsmp;
    fprintf(stderr, "nominal segment length = %d samples (%.6f sec)\n", ntake, ttake);

    /* Compute nominal number {njump} of sampling steps between seg starts: */
    int njump = (int)floor(o->every * spu + 0.5);
    if (njump == 0) { njump = 1; }
    assert(njump > 0);
    double tjump = ((double)njump)/fsmp;
    fprintf(stderr, "nominal segment stride = %d samples (%.6f sec)\n", njump, tjump);
    
    /* Compute fractional half-extent {rjoin} of splice, in sampling steps: */
    double rjoin = o->splice/2 * spu;
    assert(rjoin >= 0.0);
    double tjoin = rjoin/fsmp;
    fprintf(stderr, "splicing radius = %.6f samples (%.6f sec)\n", rjoin, tjoin);

    /* Compute extra samples {njoin} to take on internal cuts: */
    int njoin = (int)ceil(rjoin + 0.5) - 1; 
    assert(njoin >= 0);
    assert(rjoin <= njoin + 0.5);
    
    /* Compute number {nsegs} of segments to take. */
    int nsegs = (ns < ntake ? 0 : (ns - ntake) / njump + 1);
    if (nsegs == 0) { fprintf(stderr, "input file is too short"); }
    fprintf(stderr, "extracting %d segments\n", nsegs);
    
    /* Compute number {nspan} of samples spanned by all segs in input stream: */
    int nspan = (nsegs == 0 ? 0 : ntake + njump * (nsegs - 1));
    assert(ns >= nspan);
    
    /* Compute index {iskip} of first sample of first segment in input: */
    int iskip = (ns - nspan) / 2;
    
    /* Compute index {oskip} of that sample in output: */
    int oskip = 0;
    
    /* Compute number of samples in output file: */
    int nosmp = nsegs * ntake;
    fprintf(stderr, "total samples in output file = %d\n", nosmp);
    
    /* Allocate and and initialize a sound clip for the output signal: */
    sound_t so = jsa_allocate_sound(nc, nosmp);
    so.fsmp = fsmp;
    so.ns = nosmp;
    
    /* Next we extract the segments from the input file and
      concatenate them. Instead of reading the whole input file into
      memory, we read only those chunks that contain the requested
      segments, plus any samples needed for the smooth joins. */
    
    /* Allocate a sound clip structure large enough for the largest input chunk: */
    int nfull = ntake + 2*njoin;
    sound_t si = jsa_allocate_sound(nc, nfull);
    si.fsmp = fsmp;
    si.ns = nfull;

    /* Position and size of the current chunk: */
    int ichunk = -1; /* Index of sample {si[0]} in input file. */
    int nchunk = -1; /* Number of samples used in {si}. */

    /* Loop on segments/chunks: */
    int iseg; /* Index of segment. */
    int nread = 0; /* Number of samples per channel read from the input file. */
    for (iseg = 0; iseg < nsegs; iseg++)
      { 
        /* The previous input chunk, if any, is stored in {si[0..nchunk]},
          and sample {si[0]} has index {ichunk} in the input file. */
        
        /* Compute smooth join radii on each side: */
        double rlo = (iseg == 0 ? 0.0 : rjoin);
        double rhi = (iseg == nsegs - 1 ? 0.0 : rjoin);
        
        /* Compute number of extra samples to take on each side: */
        int nlo = (int)ceil(rlo + 0.5) - 1; 
        int nhi = (int)ceil(rhi + 0.5) - 1; 

        /* Save data of previous chunk: */
        int ichunk_old = ichunk;
        int nchunk_old = nchunk;
        
        /* Compute index {ichunk} of first sample of chunk in input file: */
        int ichunk = iskip + iseg*njump - nlo;

        /* Compute number {nchunk} of segment samples to read: */
        int nchunk = nlo + ntake + nhi;
        
        /* Read the relevant chunk from the input file: */
        if (ichunk >= nread)
          { /* Chunks do not overlap. */
            int nskip = ichunk - nread;
            fprintf(stderr, "skipping %d samples [%d..%d]\n", nskip, nread, nread + nskip - 1);
            jsa_skip_au_file_samples(stdin, &h, nc*nskip); 
            nread += nskip;
            fprintf(stderr, "reading %d samples [%d..%d]\n", nchunk, nread, nread + nchunk - 1);
            jsa_read_au_file_samples(stdin, &h, &si, 0, nchunk);
            nread += nchunk;
          }
        else
          { /* This chunk overlaps the previous one. */
            assert(nchunk_old >= 0);
            assert(ichunk_old >= 0);
            assert(ichunk >= nread - nchunk_old);
            /* Shift down the overlapped portion in {si}: */
            int nover = nread - ichunk;
            int nskip = nchunk_old - nover;
            int k;
            fprintf(stderr, "reusing %d samples [%d..%d]\n", nover, nread - nover, nread - 1);
            for (k = 0; k < nover; k++)
              { int c; 
                for (c = 0; c < nc; c++)
                  { si.sv[c][k] = si.sv[c][nskip + k]; }
              }
            /* Read the new portion: */
            int nrest = nchunk - nover;
            fprintf(stderr, "reading %d samples [%d..%d]\n", nrest, nread, nread + nrest - 1);
            jsa_read_au_file_samples(stdin, &h, &si, nover, nrest);
            nread += nrest;
          }
        
        /* Compute index {opos} of that sample in output clip: */
        int opos = oskip + iseg*ntake;

        /* Cut and paste segment (with smooth join if applicable): */
        cut_and_paste_segment(&si, nlo,  ntake, rlo, rhi,  &so, opos);
      }
      
    /* Output sound: */
    write_signal(stdout, &so);

    return 0;
  }

void write_signal(FILE *wr, sound_t *s)
  { fprintf(stderr, "Writing output file");
    fprintf(stderr, " (Sun \".au\" format)\n");
    jsa_write_au_file(wr, s);
  }
  
double samples_per_unit(time_unit_t unit, double fsmp)
  {
    if (unit == TUN_SECOND)
      { return fsmp; }
    else if (unit == TUN_SAMPLE)
      { return 1.0; }
    else
      { assert(FALSE); }
  }

void cut_and_paste_segment
  ( sound_t *si, 
    int ipos, 
    int ntake, 
    double rlo, 
    double rhi,
    sound_t *so, 
    int opos
  )
  {
    /* Number of channels: */
    int nc = (si->nc < so->nc ? si->nc : so->nc);
    
    /* Number of extra splicing samples before segment: */
    int nlo = (int)ceil(rlo + 0.5) - 1;
    assert(rlo <= nlo + 1);
    assert(ipos - nlo >= 0); 
    assert(opos - nlo >= 0); 
    
    /* Number of extra splicing samples after segment: */
    int nhi = (int)ceil(rhi + 0.5) - 1;
    assert(rhi <= nhi + 1);
    assert(ipos + ntake + nhi <= si->ns); 
    assert(opos + ntake + nhi <= so->ns);

    /* Loop on channels: */
    int c;
    for (c = 0; c < nc; c++)
      { int k;
        for (k = -nlo; k < ntake + nhi; k++) 
          { /* Extract sample from input signal: */
            double sk = si->sv[c][ipos+k];
            /* Compute clip mask {wk}: */
            double wlo = smooth_step(k + 0.5, rlo);
            double whi = smooth_step(ntake - k - 0.5, rhi);
            /* Paste onto the output signal: */
            so->sv[c][opos+k] += wlo*whi*sk;
          }
      }
  }

double smooth_step(double x, double r)
  { if (x <= -r) 
      { return 0.0; }
    else if (x >= +r)
      { return 1.0; }
    else
      { /* Hann half-window: */
        return (1 + sin(M_PI*(x/r)/2))/2;
      }
  }

options_t* get_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    argparser_get_keyword(pp, "-take");
    o->take = argparser_get_next_double(pp, 1.0e-100, DBL_MAX);
    
    argparser_get_keyword(pp, "-every");
    o->every = argparser_get_next_double(pp, o->take, DBL_MAX);
    
    if (argparser_keyword_present(pp, "-splice"))
      { o->splice = argparser_get_next_double(pp, 0.0, DBL_MAX); }
    else 
      { o->splice = 0.0; }
      
    if (argparser_keyword_present(pp, "-unit"))
      { 
        if (argparser_keyword_present_next(pp, "second"))
          { o->unit = TUN_SECOND; }
        else if (argparser_keyword_present_next(pp, "sample"))
          { o->unit = TUN_SAMPLE; }
        else
          { argparser_error(pp, "time unit is missing or invalid"); }
      }
    else 
      { o->unit = TUN_SECOND; }
    
    argparser_finish(pp);
    return o;
  }
