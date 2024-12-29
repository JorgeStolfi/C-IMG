#define PROG_NAME "aufilter"
#define PROG_DESC "frequency domain filtering of a audio file"
#define PROG_VERS "1.0"

/* Copyright © 2006 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2024-12-21 14:01:50 by stolfi */

#define PROG_HELP \
  PROG_NAME "\\\n" \
  "    -window {WINDOW_SIZE}  [ -overlap {WINDOW_OVERLAP} ] \\\n" \
  "    [ -noiseInput ] \\\n" \
  "    [ -fmin {FREQ_MIN} ] \\\n" \
  "    [ -fmax {FREQ_MAX} ] \\\n" \
  "    [ -kill {MIN_GAIN} {FREQ_CTR} {FREQ_DEV} ]... \\\n" \
  "    [ -chaff ] \\\n" \
  "    [ -debugTime {DEBUG_TIME} ] \\\n" \
  "    " argparser_help_info_HELP "\\\n" \
  "    < {INPUT_FILE} \\\n" \
  "    > {OUTPUT_FILE} \\\n"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads an audio file and applies a" \
  " frequency-domain (linear, time-invariant) filter to it.\n" \
  "\n" \
  "  The input signal is read from standard input and is converted from" \
  " its external format to a sequence of double-precision floating-point" \
  " samples.  This sequence is processed in overlapping batches (frames)" \
  " of {WINDOW_SIZE} consecutive samples, where successive batches are" \
  " displaced by {WINDOW_SIZE}/{WINDOW_OVERLAP} samples.  The program" \
  " multiplies each" \
  " frame by a smooth partition-of-unity windowing function, computes" \
  " its discrete Fourier transform (DFT), applies the specified" \
  " filter, computes the inverse transform, and" \
  " adds it to the output signal, which is written to the output file.\n" \
  "\n" \
  "  The filter may include a bandpass filter and zero or more" \
  " band-kill filters.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -window {WINDOW_SIZE} \n" \
  "    Specifies the number {WINDOW_SIZE} of samples per frame.  Values" \
  " of {WINDOW_SIZE} between 1024 and 4096 are usually adequate for human" \
  " speech and music (at 44100 Hz sampling rate).  This parameter is required.\n" \
  "\n" \
  "  -overlap {WINDOW_OVERLAP}\n" \
  "    Specifies the number of consecutive frames that overlap." \
  " This value must be positive and even.  The default is 8.\n" \
  "\n" \
  "  -fmin {FREQ_MIN}\n" \
  "    Specifies the cutoff frequency of the highpass" \
  " filter.  The default is 0 (meaning no filtering is done).\n" \
  "\n" \
  "  -fmax {FREQ_MAX}\n" \
  "    Specifies the cutoff frequency of the lowpass" \
  " filter.  The default is infinity (meaning no filtering is done).\n" \
  "\n" \
  "  -kill {MIN_GAIN} {FREQ_CTR} {FREQ_RAD}\n" \
  "    Pipes the signal through a bandkill filter with complemented Gaussian" \
  " profile (in log frequency space), center frequency {FREQ_CTR}" \
  " and root-mean-square band half-width {FREQ_DEV}.  The minimum gain" \
  " will be {MIN_GAIN}, at {FREQ_CTR}.  The default" \
  " is no bandpass filter.  This parameter may be repeated; each" \
  " instance specifies an additional filter stage.\n" \
  "\n" \
  "  -chaff\n" \
  "    If present, this otional parameter specifies that the output must be the noise that would" \
  " be removed, instead of the cleaned-up signal. \n" \
  "\n" \
  "  -noiseInput\n" \
  "    If present, this otional parameter causes the input signal to be replaced" \
  " by a synthetic white noise signal of the same length and sampling frequency.\n" \
  "\n" \
  "  -debugTime {DEBUG_TIME}\n" \
  "    If present, this otional parameter asks for debugging" \
  " output to be produced for the frame closest to the specified" \
  " time coordinate (in seconds). \n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  spegram(1), ausample(1), auclean(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created fev/2023 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP), by splitting out the linear" \
  " filtering part of {auclean.c}."

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
#include <vec.h>
#include <nget.h>
#include <fget.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <rn.h>
#include <argparser.h>
#include <jsrandom.h>

#include <jsaudio.h>
#include <jsaudio_au.h>

/* DATA TYPES */

/* Command line options: */
typedef struct options_t
  { /* Input parameters: */
    bool_t noiseInput;    /* TRUE replaces the input signal by white noise. */
    /* Fourier analysis parameters: */
    uint32_t wsize;        /* Window width. */
    uint32_t overlap;      /* Num of consecutive overlapping frames. */
    /* Parameters of the bandpass filter: */
    double fmin;          /* Cutoff frequency of the highpass filter. */
    double fmax;          /* Cutoff frequency of the lowpass filter. */
    double_vec_t gctr;    /* Minimum gains of the bandkill filters. */
    double_vec_t fctr;    /* Center frequencies of the bandkill filters. */
    double_vec_t fdev;    /* Nominal half-widths of the bandkill filters. */
    /* Output parameters: */
    bool_t chaff;         /* TRUE discards the wheat and keeps the chaff. */
    /* Debugging parameters: */
    double debugTime;     /* Time coordinate of frame to debug, or negative if none. */
  } options_t;

/* PROTOTYPES */

int32_t main(int32_t argc, char **argv);
  /* Main program. */

jsaudio_t read_signal(char *fname, char *sname);
  /* Reads a sound clip from the file named {fname} (or from {stdin},
    if {fname} is "stdin"). The file must be in the Sun ".au" audio file format.
    The {sname} is the nature of the signal, e.g. "input" or "noise sample". */

void write_signal(FILE *wr, jsaudio_t *s);
  /* Writes a sound clip to stream {wr}, in the Sun ".au" audio file format. */

void combine_fixed_filters
  ( uint32_t nw,           /* Number of Fourier components. */
    double fstep,         /* Frequency step (Hz). */
    double fmin,          /* Cutoff frequency of the highpass filter (Hz). */
    double fmax,          /* Cutoff frequency of the lowpass filter (Hz). */
    double_vec_t *gctr,   /* Central gains of the bandkill filters. */
    double_vec_t *fctr,   /* Center frequencies of the bandkill filters (Hz). */
    double_vec_t *fdev,   /* Nominal half-widths of the bandkill filters (Hz). */
    double gain[]         /* OUT: Transfer function of combined filters. */
  );
  /* Computes the real transfer function {gain[0..nw-1]} resulting
    from a highpass filter with cutoff frequency {fmin}, a lowpass filer
    with cutoff frequency {fmax}, and zero or more bandkill filters
    with center gains {gctr[j]}, center frequencies {fctr[j]}, and
    half-widths {fdev[j]}. Assumes step {fstep} between consecutive
    frequencies. */

void compose_bandpass_filter
  ( uint32_t nw, 
    double fstep, 
    double fmin, 
    double fmax, 
    double gain[]
  );
  /* Multiplies the real transfer function {gain[0..nw-1]} by a
    bandpass filter with cutoff frequencies {fmin} and {fmax}. Assumes
    step {fstep} between consecutive frequencies. */

double highpass_filter_gain(double freq, double fmin, double fdev);
  /* Computes the gain at frequency {freq} of a highpass filter with cutoff frequency
    {fmin} (both in Hz). The filter has smooth roll-off with nominal half-width
    {fdev} (Hz). */

void compose_bandkill_filter
  ( uint32_t nw,
    double fstep,
    double gctr, 
    double fctr, 
    double fdev, 
    double gain[]
  );
  /* Multiplies the real transfer function {gain[0..nw-1]} by the
    transfer functions of zero or more bandkill filters with center
    gains {gctr[j]}, center frequencies {fctr[j]}, and half-widths
    {fdev[j]}. Assumes step {fstep} between consecutive
    frequencies. */

void process_frame
  ( jsaudio_t *si,          /* Input sound signal. */
    uint32_t skip,         /* Index of first sample of frame yo be processed. */
    uint32_t nw,           /* Samples per frame (window). */
    uint32_t no,           /* Number of successive overlapping frames. */
    double gain[],        /* Transfer function of fixed filters. */
    fftw_complex in[],    /* WORK: FFT input vector. */
    fftw_complex out[],   /* WORK: FFT output vector. */
    fftw_plan *plan_dir,  /* Precomputed direct FFTW parameters ({in ==> out}). */
    fftw_plan *plan_inv,  /* Precomputed inverse FFTW parameters ({out ==> in}). */
    jsaudio_t *so,          /* Output sound signal. */
    bool_t debug          /* Generate diagnostic output. */
  );
  /* Processes one frame of the input signal {si}, consisting of {nw}
    samples starting at sample {skip}, and pastes the filtered result
    onto the output signal {so}. Assumes that the frame overlaps the
    next {no-1} frames and has been apodized by the proper windowing
    function.
    
    The processing consists of applying a fixed filter with the
    transfer function {gain}. The vectors {in},
    {out}, and {totgain} must be allocated by the caller. */
    
void apply_hann_window_to_signal(fftw_complex *sg, uint32_t nw, uint32_t no);
  /* Multiplies the signal {sg[0..nw-1]} by a Hann smoothing window,
    scaled so that the windows of all frames, assuming {no} overlap factor,
    add up to the unit constant function. */
    
void show_filter(uint32_t nw, double gain[]);
  /* Writes to {stderr} the transfer function {gain[0..nw]} of some filter. */

options_t* get_options(int32_t argc, char **argv);
  /* Parses the command-line arguments. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    /* Read the input signal: */
    jsaudio_t si = read_signal("stdin", "input");
    
    uint32_t nw = o->wsize;    /* Number of samples in window. */
    fprintf(stderr, "window size = %d samples (%.6f sec)\n", nw, ((double)nw)/si.fsmp);
    uint32_t nfrq = nw/2 + 1;   /* Number of distinct frequencies, per Nyquist. */
    fprintf(stderr, "spectrum has %d frequencies {0..%d}\n", nfrq, nfrq-1);
    uint32_t no = o->overlap;  /* Number of consecutive samples that overlap. */
    fprintf(stderr, "frame overlap count = %d\n", no);
    assert(no % 2 == 0);
    assert(nw % no == 0);
    uint32_t stride = nw/no; /* Number of time steps between spectrograms. */
    uint32_t nfrm = (si.ns - nw)/stride + 1;  /* Number of frames in input signal. */
    fprintf(stderr, "number of frames = %d\n", nfrm);
    fprintf(stderr, "frame spacing = %d samples (%.6f sec)\n", stride, ((double)stride)/si.fsmp);
    double fstep = si.fsmp/nw;  /* Difference between consecutive frequencies. */
    
    /* If requested, substitute white noise for the input: */
    if (o->noiseInput)
      { /* Fill frame with white noise: */
        for (uint32_t c = 0;  c < si.nc; c++)
          { for (uint32_t ismp = 0;  ismp < si.ns; ismp++) 
             { si.sv[c][ismp] = 2*drandom() - 1; }
          }
      }
    
    /* Allocate working storage for Fast Fourier Transform (FFT): */
    fftw_complex in[nw];  /* FFTW input vector. */
    fftw_complex out[nw]; /* FFTW output vector. */

    /* Precompute parameters for the FFT: */
    fftw_plan plan_dir = fftw_plan_dft_1d((int32_t)nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_inv = fftw_plan_dft_1d((int32_t)nw, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Allocate transfer function (TFs) of fixed filters: */
    double fixgain[nw]; /* TF of all fixed filters (precomputed). */

    /* Precompute TF of fixed filters: */
    combine_fixed_filters
      ( nw, fstep,
        o->fmin, o->fmax, 
        &(o->gctr), &(o->fctr), &(o->fdev), 
        fixgain
      );

    /* If user wants the removed chaff, complement the gain: */
    if (o->chaff) 
      { for (uint32_t k = 0;  k < nw; k++) { fixgain[k] = 1.0 - fixgain[k]; }}

    /* Allocate and and initialize the output signal: */
    jsaudio_t so = jsaudio_allocate_sound(si.nc, si.ns);
    so.fsmp = si.fsmp;
    so.ns = si.ns;
    
    /* Compute frames and output them: */
    double totTime = si.ns/((double)si.fsmp);
    int32_t debug_ifrm = (int32_t)floor(nfrm*(o->debugTime/totTime));
    for (uint32_t ifrm = 0;  ifrm < nfrm; ifrm++)
      { /* Compute index {skip} of first sample in frame: */
        uint32_t skip = ifrm*stride;
        /* Compute the power spectrum: */
        bool_t debug = (ifrm == debug_ifrm);
        process_frame
          ( &si, skip, 
            nw, no, 
            fixgain, 
            in, out, &plan_dir, &plan_inv,
            &so,
            debug
          );
      }
      
    /* Output sound: */
    write_signal(stdout, &so);

    return 0;
  }
  
jsaudio_t read_signal(char *fname, char *sname)
  { fprintf(stderr, "Reading %s file from %s...\n", sname, fname);
    FILE *rd = (strcmp(fname, "stdin") == 0 ? stdin : open_read(fname, TRUE));
    jsaudio_t si = jsaudio_au_read_file(rd);
    fprintf(stderr, "input channels = %d\n", si.nc);
    fprintf(stderr, "input samples %d\n", si.ns);
    fprintf(stderr, "sampling frequency = %.4f  Hz\n", si.fsmp);
    double tstep = 1.0/si.fsmp;
    fprintf(stderr, "sampling step = %.6f sec\n", tstep);
    double totTime = si.ns*tstep;
    fprintf(stderr, "total duration = %.6f sec\n", totTime);
    if (rd != stdin) { fclose(rd); }
    return si;
  }

void write_signal(FILE *wr, jsaudio_t *s)
  { fprintf(stderr, "Writing output file");
    fprintf(stderr, " (Sun \".au\" format)\n");
    jsaudio_au_write_file(wr, s);
  }
  
void combine_fixed_filters
  ( uint32_t nw,           /* Number of Fourier components. */
    double fstep,         /* Frequency step (Hz). */
    double fmin,          /* Cutoff frequency of the highpass filter (Hz). */
    double fmax,          /* Cutoff frequency of the lowpass filter (Hz). */
    double_vec_t *gctr,   /* Central gains of the bandkill filters. */
    double_vec_t *fctr,   /* Center frequencies of the bandkill filters (Hz). */
    double_vec_t *fdev,   /* Nominal half-widths of the bandkill filters (Hz). */
    double gain[]          /* OUT: Transfer function of combined filters. */
  )
  {
    /* Initialize the filter's transfer function: */
    for (uint32_t k = 0;  k < nw; k++) { gain[k] = 1.0; }

    /* Compose the highpass and lowpass filters: */
    if ((fmin > 0.0) || (fmax < INFINITY))
      { compose_bandpass_filter(nw, fstep, fmin, fmax, gain); }

    /* Add effect of bandkill filters: */
    assert(fctr->ne == fdev->ne);
    uint32_t nk = fctr->ne;
    for (uint32_t j = 0;  j < nk; j++)
      { compose_bandkill_filter
          ( nw, fstep,
            gctr->e[j], fctr->e[j], fdev->e[j], 
            gain
          );
      }
    for (uint32_t k = 0;  k < nw; k++) 
      { assert(isfinite(gain[k]));
        assert(gain[k] >= 0);
      }
  }

void compose_bandpass_filter
  ( uint32_t nw, 
    double fstep, 
    double fmin, 
    double fmax, 
    double gain[]
  )
  {
    /* Number of frequencies: */
    uint32_t nfrq = nw/2 + 1;
    
    for (uint32_t fp = 0;  fp < nfrq; fp++)
      { /* Get the freq {fm} equivalent to freq {-fp} modulo {nw}: */
        uint32_t fm = (nw - fp) % nw; 
        /* Absolute frequency in Hertz: */
        double freq = fp * fstep;
        /* Highpass and lowpass gains: */
        double logain = highpass_filter_gain(freq, fmin, 1.0*fstep);
        double higain = 1.0 - highpass_filter_gain(freq, fmax, 3.0*fstep);
        /* Apply to TF: */
        gain[fp] *= logain*higain;
        if (fm != fp) { gain[fm] *= logain*higain; }
      }
  }
  
double highpass_filter_gain(double freq, double fmin, double fdev)
  {
    freq = fabs(freq);
    /* Check for trivial filters: */
    if (fmin == INFINITY) { return 0.0; }
    if (fmin == 0.0) { return 1.0; }
    /* Filter is not trivial; check for extreme freqs: */
    if (freq == 0.0) { return 0.0; }
    if (freq == INFINITY) { return 1.0; }
    /* The following code should be safe. */
    double Lf = log(freq);
    double Lm = log(fmin);
    double rdev = fdev/fmin;
    if (Lf < Lm - 5.0*rdev) { return 0.0; }
    if (Lf > Lm + 5.0*rdev) { return 1.0; }
    double xi = (Lf - Lm)/rdev;
    return (1 + erf(xi))/2;
  }
  
void compose_bandkill_filter
  ( uint32_t nw,
    double fstep,
    double gctr, 
    double fctr, 
    double fdev, 
    double gain[]
  )
  {
    /* Number of frequencies: */
    uint32_t nfrq = nw/2 + 1;
    
    /* Cleanup the spectrum: */
    for (uint32_t fp = 0;  fp < nfrq; fp++)
      { /* Get the freq {g} equivalent to freq {-f}. */
        uint32_t fm = (nw - fp) % nw;  /* Frequency {-fp} modulo {nw}. */

        if (fp > 0) 
          { /* Compute the absolute frequency {freq} in Hertz: */
            double freq = fp*fstep;

            /* Go to log-frequency space: */
            double logfreq = log(freq);
            /* Get log center and log half-width of bandkill filter {j}: */
            double logfctr = log(fctr);
            double rdev = fdev/fctr;
            double off = logfreq - logfctr;
            if (fabs(off) < 5.0*rdev)
              { /* Worth computing the Gaussian: */
                double x = off/rdev;
                double w = 1.0 - (1.0 - gctr)*exp(-x*x/2);
                gain[fp] *= w;
                if (fm != fp) { gain[fm] *= w; }
              }
          }
      }
  }

void process_frame
  ( jsaudio_t *si,          /* Input sound signal. */
    uint32_t skip,         /* Index of first sample of frame yo be processed. */
    uint32_t nw,           /* Samples per frame (window). */
    uint32_t no,           /* Number of successive overlapping frames. */
    double fixgain[],     /* Transfer function of fixed filters. */
    fftw_complex in[],    /* WORK: FFT input vector. */
    fftw_complex out[],   /* WORK: FFT output vector. */
    fftw_plan *plan_dir,  /* Precomputed direct FFTW parameters ({in ==> out}). */
    fftw_plan *plan_inv,  /* Precomputed inverse FFTW parameters ({out ==> in}). */
    jsaudio_t *so,          /* Output sound signal. */
    bool_t debug          /* Generate diagnostic output. */
  )
  { 
    assert(si->fsmp == so->fsmp);
    
    /* Number of channels: */
    uint32_t nc = (si->nc < so->nc ? si->nc : so->nc);
    
    /* Seed for the noise generator: */
    srandom(519257027);
    
    /* Loop on channels: */
    for (uint32_t c = 0;  c < nc; c++)
      { /* Extract frame from input signal: */
        assert(skip + nw <= si->ns); 
        for (uint32_t k = 0;  k < nw; k++) 
          { in[k][0] = si->sv[c][skip+k]; in[k][1] = 0.0; }
        
        /* Apply the Hann unit-partition window function: */
        apply_hann_window_to_signal(in, nw, no); 
          
        /* Compute direct FFT: */
        fftw_execute(*plan_dir);
       
        /* if (skip == 0) { show_filter(nw, fixgain); } */
        
        /* Apply filters to the frame spectrum: */
        for (uint32_t k = 0;  k < nw; k++) 
          { out[k][0] *= fixgain[k]; out[k][1] *= fixgain[k]; }
        
        /* Compute inverse FFT: */
        fftw_execute(*plan_inv);

        /* Splat the frame onto the output signal: */
        assert(skip + nw <= so->ns); 
        for (uint32_t k = 0;  k < nw; k++) 
          { in[k][0] /= nw;
            in[k][1] /= nw;
            so->sv[c][skip+k] += in[k][0]; 
            /* Result must be real: */
            if(fabs(in[k][1]) > 1.0e-6)
              { fprintf(stderr, "%4d %16.10f %16.10f\n", k, in[k][0], in[k][1]); } 
          }
      }
  }

void show_filter(uint32_t nw, double gain[])
  { for (uint32_t k = 0;  k < nw; k++) 
      { fprintf(stderr, "gain[%03d] = %16.12f\n", k, gain[k]); }
  }

void apply_hann_window_to_signal(fftw_complex *sg, uint32_t nw, uint32_t no)
  {
    /* Scaling factor for Hann unit-partition window function */
    assert(nw % no == 0);
    double fac = 1.0/((double)no);

    for (uint32_t k = 0;  k < nw; k++) 
      { double x = M_PI*(2*((double)k)/((double)nw) - 1.0);
        double w = fac*(1.0 + cos(x)); 
        sg[k][0] *= w; sg[k][1] *= w;
      }
  }
#define MAX_WSIZE 9192

options_t* get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    o->noiseInput = argparser_keyword_present(pp, "-noiseInput");
    
    argparser_get_keyword(pp, "-window");
    o->wsize = (uint32_t)argparser_get_next_int(pp, 1, MAX_WSIZE);
    
    if (argparser_keyword_present(pp, "-overlap"))
      { o->overlap = (uint32_t)argparser_get_next_int(pp, 2, MAXINT); }
    else 
      { o->overlap = 8; }
    if (o->overlap % 2 != 0) 
      { argparser_error(pp, "overlap count must be even"); }
    if (o->wsize % o->overlap != 0) 
      { argparser_error(pp, "overlap count must divide frame size"); }
      
    if (argparser_keyword_present(pp, "-fmin"))
      { o->fmin = argparser_get_next_double(pp, 0, DBL_MAX); }
    else 
      { o->fmin = 0.0; }
    
    if (argparser_keyword_present(pp, "-fmax"))
      { o->fmax = argparser_get_next_double(pp, o->fmin, DBL_MAX); }
    else 
      { o->fmax = INFINITY; }
    
    o->gctr = double_vec_new(0);
    o->fctr = double_vec_new(0);
    o->fdev = double_vec_new(0);
    uint32_t nk = 0;
    while(argparser_keyword_present(pp, "-kill"))
      { double_vec_expand(&(o->gctr), nk);
        o->gctr.e[nk] = argparser_get_next_double(pp, 0, DBL_MAX);
        double_vec_expand(&(o->fctr), nk);
        o->fctr.e[nk] = argparser_get_next_double(pp, 1.0e-100, 1.0e+100);
        double_vec_expand(&(o->fdev), nk);
        o->fdev.e[nk] = argparser_get_next_double(pp, 1.0e-100, 1.0e+100);
        nk++;
      }
    double_vec_trim(&(o->gctr), nk);
    double_vec_trim(&(o->fctr), nk);
    double_vec_trim(&(o->fdev), nk);
    
    o->chaff = argparser_keyword_present(pp, "-chaff");
    
    if (argparser_keyword_present(pp, "-debugTime"))
      { o->debugTime = argparser_get_next_double(pp, -DBL_MAX, DBL_MAX); }
    else 
      { o->debugTime = -1.0; }
    
    argparser_finish(pp);
    return o;
  }
