#define PROG_NAME "ausample"
#define PROG_DESC "extract sample clips from audio files"
#define PROG_VERS "1.0"

/* Copyright © 2006 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2024-12-21 14:01:44 by stolfi */

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
  " samples.\n" \
  "\n" \
  "  Each segment to be extracted has length {CLIP_LENGTH}, and" \
  " successive segments start {CLIP_STRIDE} apart.  The series" \
  " is centered over the length of" \
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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <values.h>
#include <assert.h>
#include <math.h>

#include <fftw3.h>

#include <bool.h>
#include <jsprintf.h>
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
    TUN_SAMPLE   /* Unit is the sampling step. */
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

int32_t main(int32_t argc, char **argv);
  /* Main program. */

void write_signal(FILE *wr, jsaudio_t *s);
  /* Writes a sound clip to stream {wr}, in the Sun ".au" audio file format. */

double samples_per_unit(time_unit_t unit, double fsmp);
  /* Computes the number of sampling steps coresponding to the {unit},
    assuming that the sampling frequency is {fsmp} (in Hz). */

void paste_segment
  ( jsaudio_t *si,
    uint32_t ni,
    uint32_t nlo, 
    jsaudio_t *so, 
    uint32_t *no_P
  );
  /* Pastes {ni} samples {si[0..ni-1]} it onto {so}. 
    Assumes that, on input, {so} has {no}
    samples, where {no} is {*no_P}.  The first {nlo} samples of 
    {si} will be mixed with the last {nlo} samples
    of {so}, and the rest will be stored into {so[no..no_new-1]}
    where {no_new = no+ni-nlo}.  Updates {*no_P} with {no_new}.
    
    The number {nlo} must not exceed {no}. */
 
double splice_weight(uint32_t k, uint32_t n);
  /* Computes the weight of sample {k} of the input 
    in a splice that involves the first {n} samples.
    Requires {k} in {0..n-1}. */

options_t* get_options(int32_t argc, char **argv);
  /* Parses the command-line arguments. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    /* Read the input file header: */
    jsaudio_au_file_header_t h = jsaudio_au_read_file_header(stdin);
    
    /* Compute number of channels {nc} and samples per channel {ns}: */
    uint32_t nc = h.channels;
    uint32_t bps = jsaudio_au_file_bytes_per_sample(h.encoding);
    uint32_t ns = h.data_size / nc / bps;
    assert(h.data_size == ns * nc * bps); 
    fprintf(stderr, "input channels = %d\n", nc);
    fprintf(stderr, "input samples %d\n", ns);
    
    /* Get sampling frequency: */
    double fsmp = (double)h.sample_rate;
    fprintf(stderr, "sampling frequency = %.4f  Hz\n", fsmp);
    
    /* Get number of samples per unit: */
    double spu = samples_per_unit(o->unit, fsmp);

    /* Compute nominal number {ntake} of samples per segment (excluding splices): */
    uint32_t ntake = (uint32_t)floor(o->take * spu + 0.5); 
    demand(ntake >= 0, "segments have zero samples");
    double ttake = ((double)ntake)/fsmp;
    fprintf(stderr, "nominal segment length = %d samples (%.6f sec)\n", ntake, ttake);

    /* Compute nominal number {njump} of sampling steps between seg starts: */
    uint32_t njump = (uint32_t)floor(o->every * spu + 0.5);
    if (njump == 0) { njump = 1; }
    double tjump = ((double)njump)/fsmp;
    fprintf(stderr, "nominal segment stride = %d samples (%.6f sec)\n", njump, tjump);
    
    /* Compute extra samples {njoin} for splicing: */
    double rjoin = o->splice/2 * spu;
    assert(rjoin >= 0.0);
    uint32_t njoin = (uint32_t)ceil(rjoin + 0.5) - 1; 
    assert(njoin >= 0);
    assert(rjoin <= njoin + 0.5);
    double tjoin = njoin/fsmp;
    fprintf(stderr, "splice samples = %d (%.6f sec)\n", njoin, tjoin);
    
    /* Compute number {nsegs} of segments to take. */
    uint32_t nsegs;
    if (ns < ntake)
      { nsegs = 0; }
    else if (ns < 2*ntake + njump)
      { nsegs = 1; }
    else
      { nsegs = (ns - ntake) / njump + 1;
        demand(njump >= ntake + 2*njoin, "segment tails overlap");
      }
        
    demand(nsegs >= 1, "input file is too short");
    fprintf(stderr, "extracting %d segments\n", nsegs);
    
    /* Compute number {nspan} of samples spanned by all segs in input stream: */
    uint32_t nspan = (nsegs == 1 ? ntake : ntake + njump*(nsegs-1));
    fprintf(stderr, "reading a total of %d samples\n", nspan);
    assert(ns >= nspan);
    
    /* Compute index {iskip} of first sample of first segment in input: */
    uint32_t iskip = (ns - nspan) / 2;
    
    /* Compute number of samples in output file: */
    uint32_t nosmp = (nsegs == 1 ? ntake : ntake*nsegs + njoin*(nsegs-1));
    fprintf(stderr, "total samples in output file = %d\n", nosmp);
    
    /* Allocate and and initialize a sound clip for the output signal: */
    jsaudio_t so = jsaudio_allocate_sound(nc, nosmp);
    so.fsmp = fsmp;
    so.ns = nosmp;
    
    /* Next we extract the segments from the input file and
      concatenate them. Instead of reading the whole input file into
      memory, we read only those chunks that contain the requested
      segments, plus any samples needed for the smooth joins. */
    
    /* Allocate a {jsaudio_t} large enough for the largest input chunk: */
    uint32_t nfull = ntake + 2*njoin;
    jsaudio_t si = jsaudio_allocate_sound(nc, nfull);
    si.fsmp = fsmp;
    si.ns = nfull;

    /* Loop on segments/chunks: */
    uint32_t nread = 0; /* Number of samples per channel read from the input file. */
    uint32_t no = 0; /* Samples in {so}. */
    for (uint32_t iseg = 0; iseg < nsegs; iseg++)
      { 
        /* The previous segments, spliced, are stored in 
          {so[0..onext-1] where {onext} is 0 if {iseg} is zero,
          else {iseg*ntake + njoin}.  Already {nrad} samples
          have been read from the input file.*/
        
        /* Compute number of extra samples to take on each side: */
        uint32_t nlo = (iseg == 0 ? 0 : njoin); 
        uint32_t nhi = (iseg == nsegs - 1 ? 0 : njoin); 

        /* Compute number {nchunk} of segment samples to read: */
        uint32_t nchunk = nlo + ntake + nhi;
        
        /* Compute index {ichunk} of first sample of next chunk in input file: */
        int32_t ichunk = (int32_t)(iskip + iseg*njump) - (int32_t)nlo;
        assert(ichunk > 0);
        
        /* Read the next segment from the input file into {si[0..nchunk-1]}: */
        assert(ichunk >= nread); /* Chunks should not overlap.*/
        uint32_t nskip = (uint32_t)ichunk - nread;
        fprintf(stderr, "skipping %d samples [%d..%d]\n", nskip, nread, nread + nskip - 1);
        jsaudio_au_skip_file_samples(stdin, &h, nc*nskip); 
        nread += nskip;
        fprintf(stderr, "reading segment %d - %d samples [%d..%d]\n", iseg, nchunk, nread, nread + nchunk - 1);
        jsaudio_au_read_file_samples(stdin, &h, &si, 0, nchunk);
        nread += nchunk;
        
        /* Cut and paste segment (with smooth join if applicable): */
        paste_segment(&si, nchunk, nlo, &so, &no);
      }
      
    /* Output sound: */
    write_signal(stdout, &so);

    return 0;
  }

void write_signal(FILE *wr, jsaudio_t *s)
  { fprintf(stderr, "Writing output file");
    fprintf(stderr, " (Sun \".au\" format)\n");
    jsaudio_au_write_file(wr, s);
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

void paste_segment
  ( jsaudio_t *si, 
    uint32_t ni, 
    uint32_t nlo, 
    jsaudio_t *so, 
    uint32_t *no_P
  )
  {
    uint32_t no = (*no_P);
    demand(nlo <= no, "not enough output samples to overlap");
    demand(nlo <= ni, "not enough input samples to overlap");
    uint32_t ko = no - nlo; /* Nex sample to splice/store in {so}. */
    fprintf(stderr, "pasting into {so[%d..%d]}\n", ko, ko + ni - 1);
    demand(no + ni - nlo <= so->ns, "not enough space in {so}");
    
    /* Number of channels: */
    uint32_t nc = (si->nc < so->nc ? si->nc : so->nc);
    
    /* Loop on input samples: */
    for (uint32_t ki = 0; ki < ni; ki++) 
      { if (ki < nlo)
          { /* Splice: */
            assert(ko < no);
            double wik = splice_weight(ki, nlo);
            double wok = 1 - wik;
            for (uint32_t c = 0; c < nc; c++) 
              { double ski = si->sv[c][ki];
                double sko = so->sv[c][ko];
                so->sv[c][ko] = wok*sko + wik*ski; 
              }
          }
        else
          { /* Append: */
            for (uint32_t c = 0; c < nc; c++)
              { so->sv[c][ko] = si->sv[c][ki];
                no++;
              }
          }
        ko++;
      }
    assert(ko == no);
    (*no_P) = no;
  }

double splice_weight(uint32_t k, uint32_t n)
  { demand(k < n, "not a splice point");
    double arg = M_PI*(k + 0.5)/n;
    return 0.5*(1 - cos(arg));
  }

options_t* get_options(int32_t argc, char **argv)
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
