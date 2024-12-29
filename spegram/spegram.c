#define PROG_NAME "spegram"
#define PROG_DESC "compute the energy spectrogram of a signal"
#define PROG_VERS "1.0"

/* Copyright © 2006 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2024-12-21 11:58:31 by stolfi */

#define PROG_HELP \
  PROG_NAME "\\\n" \
  "    -input [ ascii | au ] \\\n" \
  "    -window [ hann | rect ] {WINDOW_SIZE}  [ -stride {WINDOW_STRIDE} ] \\\n" \
  "    -output [ ascii | pgm ] [ -range {VMAX} {VMIN} ] \\\n" \
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
  "  The program reads a sampled signal (e.g. an audio file), and outputs" \
  " its energy spectrogram.\n" \
  "\n" \
  "  The input signal is read from standard input and is converted from" \
  " its external format to a sequence of double-precision floating-point" \
  " samples.  This sequence is processed in overlapping batches (frames)" \
  " of {WINDOW_SIZE} consecutive samples, where successive batches are" \
  " displaced by {WINDOW_STRIDE} samples.  The program multiplies each" \
  " frame by a user-specified windowing function, computes its discrete" \
  " Fourier transform (DFT), and writes its power spectrum to the output file.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -input [ ascii | au ]\n" \
  "    Specifies the format of the input file.  The value \"ascii\" means" \
  " that the input is a text file (see below for details).  The value" \
  " \"au\" means that the input is in Sun's AU audio format (usually" \
  " identified by the \".au\" filename extension).  Only a few variants" \
  " of the AU format are presently recognized, including 16-bit signed" \
  " integer linear PCM and 32-bit floating-point linear PCM.\n" \
  "  -window [ hann | rect ] {WINDOW_SIZE} \n" \
  "    Specifies the type of the apodizing function used to avoid Gibbs" \
  " ringing, and the number {WINDOW_SIZE} of samples per frame.  Presently" \
  " recognizes the types \"hann\" (Hann shifted cosine windowing funcion) and" \
  " \"rect\" (rectangular windowing function, i.e. no apodizing).  Values" \
  " of {WINDOW_SIZE} between 1024 and 4096 are usually adequate for human" \
  " speech and music (at 44100 Hz sampling rate).\n" \
  "  -stride {WINDOW_STRIDE}\n" \
  "    Specifies the number of samples between successive" \
  " frames.  Usually this value is positive and less" \
  " than {WINDOW_SIZE}.  The default is {WINDOW_SIZE/8}.\n" \
  "\n" \
  "  -output [ ascii | pgm ]\n" \
  "    Specifies the format of the output file.  The \"ascii\" option" \
  " produces a plain text file in a format suitable for {gnuplot}.  The" \
  " \"pgm\" option produces a grayscale image in the PGM(5) format." \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "ASCII INPUT FORMAT\n" \
  "  The plain text (\"ascii\") input format begins with three lines" \
  " containing \"channels = {NUM_CHANNELS}\", \"samples = {NUM_SAMPLES}\", and" \
  " \"frequency = {SAMPLING_FREQ}\".  The sampling frequency should be" \
  " in Hertz. Then follow {NUM_SAMPLES} lines, each" \
  " containing the sampling time and {NUM_CHANNELS} sample values, one" \
  " per channel.  The values should be floating-point numbers, separated" \
  " by one or more spaces." \
  "\n" \
  "ASCII OUTPUT FORMAT\n" \
  "  The plain text (\"ascii\") output format begins with four lines" \
  " containing \"tstep = {TIME_STEP}\", \"fstep = {FREQ_STEP}\", " \
  " \"freqs = {NUM_FREQUENCIES}\", and \"times = {NUM_TIMES}\".  Then" \
  " follow {NUM_TIMES} blocks of spectrogram data, one block for each" \
  " frame in the spectrogram, separated by blank lines.  Each block" \
  " consists of {NUM_FREQUENCIES} data lines.  Each data line has" \
  " three floating-point values: the time at the center of the" \
  " frame, the frequency (a multiple of {FREQ_STEP}), and the" \
  " spectrogram's value (signal power in that frame and frequency).\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgminvert(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created oct/2006 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP), as {egrama.c}.\n" \
  "  Revised and renamed {spegram.c} on oct/2006 by Jorge Stolfi, UNICAMP."

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
#include <argparser.h>

#include <jsaudio.h>
#include <jsaudio_au.h>

/* COMMAND LINE ARGUMENTS */

/* Input file type: */
typedef enum { 
  IFM_AU,       /* Sun AU format. */
  IFM_ASCII     /* Plain ASCII format. */
} Input_Fmt_t;    

/* Window function type: */
typedef enum {
  WIN_HANN,     /* Hann window. */
  WIN_RECT      /* Rectangular window. */
} Window_Type_t; 

/* Output file type: */
typedef enum {
  OFM_PGM,     /* PGM output. */
  OFM_ASCII    /* Plain ASCII output, suitable for {gnuplot}. */
} Output_Fmt_t;   

/* Command line options: */
typedef struct options_t 
  { Input_Fmt_t ifmt;     /* Format of input file. */
    Window_Type_t wtype;  /* Window function type. */
    uint32_t wsize;            /* Window width. */
    uint32_t stride;           /* Samples to skip between successive windows. */
    double vmin;          /* Map this spectrogram value to 0.0 output. */
    double vmax;          /* Map this spectrogram value to 1.0 output. */
    Output_Fmt_t ofmt;    /* Format of output file. */
  } options_t;

/* PROTOTYPES */

jsaudio_t read_signal(FILE *rd, Input_Fmt_t ifmt);
jsaudio_t read_signal_ascii(FILE *rd);
void write_output_header(FILE *wr, Output_Fmt_t ofmt, uint32_t ne, uint32_t nf, double tstep, double fstep);
void write_frame(FILE *wr, uint32_t nfrq, double pwr[], double t, double fstep, Output_Fmt_t ofmt, double vmin, double vmax);
void write_output_trailer(FILE *wr, Output_Fmt_t ofmt);
void compute_frame
  ( jsaudio_t *s, 
    uint32_t skip, 
    Window_Type_t wtype,
    uint32_t nw,
    fftw_complex *in, 
    fftw_complex *out, 
    fftw_plan *plan,
    double pwr[]
  );
void update_egram_range(uint32_t nfrq, double pwr[], double *emin, double *emax);
double safelog(double v, double vmin);
options_t* get_options(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    jsaudio_t si = read_signal(stdin, o->ifmt);
    fprintf(stderr, "input channels = %d\n", si.nc);
    fprintf(stderr, "input samples %d\n", si.ns);
    fprintf(stderr, "sampling frequency = %.4f  Hz\n", si.fsmp);
    
    /* Sampling step: */
    double tstep = 1.0/si.fsmp;
    fprintf(stderr, "sampling step = %.6f sec\n", tstep);

    /* Get number {nw} of time steps in window: */
    uint32_t nw = o->wsize;
    fprintf(stderr, "window size = %d samples (%.6f sec)\n", nw, ((double)nw)/si.fsmp);
    
    /* Get number {stride} os samples between spectrograms: */
    uint32_t stride = o->stride;
    fprintf(stderr, "frame spacing = %d samples (%.6f sec)\n", stride, ((double)stride)/si.fsmp);
    
    /* Get number {nf} of frames in spectrogram: */
    uint32_t nf = (si.ns - nw)/stride + 1; 
    fprintf(stderr, "number of frames = %d\n", nf);
    
    /* Get number {nfrq} of values (freqs) in each frame of the spectrogram: */
    uint32_t nfrq = nw/2 + 1;
    fprintf(stderr, "number of distinct frequencies = %d\n", nfrq);
    
    /* Frequency step: */
    double fstep = si.fsmp/nw;
    fprintf(stderr, "frequency step = %.3f Hz\n", fstep);
    
    /* Allocate working storage for Fast Fourier Transform (FFT): */
    fftw_complex *in, *out;
    in = fftw_malloc(nw*sizeof(fftw_complex));     
    out = fftw_malloc(nw*sizeof(fftw_complex));    

    /* Precompute parameters for the FFT: */
    fftw_plan plan;
    plan = fftw_plan_dft_1d((int32_t)nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    /* Output header: */
    write_output_header(stdout, o->ofmt, nfrq, nf, tstep, fstep);
    
    /* Compute frames and output them: */
    double pwr[nfrq]; /* Power spectrum of slice. */
    double emin = +INFINITY;
    double emax = -INFINITY;
    for (uint32_t ifrm = 0;  ifrm < nf; ifrm++)
      { 
        /* Compute index {skip} of first sample in frame: */
        uint32_t skip = ifrm*stride;
        /* Compute time at center of window: */
        double t = (skip + 0.5*nw)/si.fsmp; 
        /* Compute the power spectrum: */
        compute_frame(&si, skip, o->wtype, nw, in, out, &plan, pwr);
        update_egram_range(nfrq, pwr, &emin, &emax);
        /* Write it out: */
        write_frame(stdout, nfrq, pwr, t, fstep, o->ofmt, o->vmin, o->vmax);
      }
      
    /* Show egram range: */
    fprintf(stderr, "spectrogram value range = [ %14.6e _ %14.6e ]\n", emin, emax); 

    /* Output header: */
    write_output_trailer(stdout, o->ofmt);

    fprintf(stderr, "done.\n");
    return 0;
  }
  
jsaudio_t read_signal(FILE *rd, Input_Fmt_t ifmt)
  { fprintf(stderr, "Reading input file");
    if (ifmt == IFM_AU)
      { fprintf(stderr, " (Sun \".au\" format)\n");
        return jsaudio_au_read_file(rd);
      }
    else if (ifmt == IFM_ASCII)
      { fprintf(stderr, " (plain ASCII format)\n");
        return read_signal_ascii(rd);
      }
    else
      { fprintf(stderr, "\n");
        fprintf(stderr, "unspecified format, aborted\n");
        exit(1);
      }
  }

jsaudio_t read_signal_ascii(FILE *rd)
  {
    uint32_t nc = nget_uint32(rd, "channels", 10); fget_eol(rd);
    uint32_t ns = nget_uint32(rd, "samples", 10);  fget_eol(rd);
    double fsmp = nget_double(rd, "frequency");  fget_eol(rd);
    jsaudio_t si = jsaudio_allocate_sound(nc, ns); 

    int32_t r = getc(rd);
    uint32_t k = 0; /* Sample index */
    while (TRUE)
      { /* Skip spaces and TABs: */
        while ((r == '\011') || (r == ' ')) { r = getc(rd); }
        /* If end-of-file, we should be done: */
        if (r == EOF) { break; }
        if (r == '#')
          { /* Comment line, ignore: */
            while ((r != EOF) && (r != '\n')) { r = getc(rd); }
          }
        else
          { /* Data line, read time {t} and {nc} channel samples: */
            ungetc(r, rd);
            (void)fget_double(rd); /* Sampling time. */
            for (uint32_t c = 0;  c < nc; c++) { si.sv[c][k] = fget_double(rd); }
            k++;
            fget_eol(rd);
          }
        r = getc(rd);
      }
    assert(k == ns);
    si.ns = ns;
    si.fsmp = fsmp;
    return si;
  }       

/* Max pixel value for PGM output: */
#define PGM_MAXVAL 65535

void write_output_header(FILE *wr, Output_Fmt_t ofmt, uint32_t nfrq, uint32_t nf, double tstep, double fstep)
  { 
    if (ofmt == OFM_PGM)
      { 
        fprintf(wr, "P2\n"); /* PGM ASCII format. */
        fprintf(wr, "%d %d\n", nfrq, nf); /* Width and height of image. */
        fprintf(wr, "%d\n", PGM_MAXVAL);
      }
    else if (ofmt == OFM_ASCII)
      {
        fprintf(wr, "# tstep = %24.16e\n", tstep);
        fprintf(wr, "# fstep = %24.16e\n", fstep);
        fprintf(wr, "# freqs = %d\n", nfrq);
        fprintf(wr, "# times = %d\n", nf);
      }
    else
      { fprintf(stderr, "unspecified format, aborted\n");
        exit(1);
      }
  }

double safelog(double v, double vmin)
  { 
    return log(hypot(v, vmin));
  }

void write_frame(FILE *wr, uint32_t nfrq, double pwr[], double t, double fstep, Output_Fmt_t ofmt, double vmin, double vmax)
  { 
    if (ofmt == OFM_PGM)
      { 
        /* Make sure that {vmin} is positive, to avoid {log} errors: */
        if (vmin <= 0) { vmin = 1.0e-100; }
        /* Ensure {vmax > vmin}, to avoid divide by zero: */
        if (vmax <= vmin) { vmax = (1 + 1.0e-14)*vmin; }
        double umin = log(vmin);
        double umax = log(vmax);
        for (uint32_t k = 0;  k < nfrq; k++)
          { if (k > 0) { fprintf(wr, "%c", (k % 10 == 0 ? '\n' : ' ')); }
            double pk = fmin(vmax, fmax(vmin, pwr[k]));
            double u = log(pk);
            double v = (u - umin)/(umax - umin); 
            assert(isfinite(v));
            if (v < 0.0) { v = 0.0; }
            if (v > 1.0) { v = 1.0; }
            uint32_t vi = (uint32_t)(v*PGM_MAXVAL + 0.5);
            assert((vi >= 0) && (vi <= PGM_MAXVAL));
            fprintf(wr, "%d", vi); 
          }
        fprintf(wr, "\n");  
      }
    else if (ofmt == OFM_ASCII)
      { for (uint32_t k = 0;  k < nfrq; k++)
          { double v = (pwr[k] - vmin)/(vmax - vmin);
            fprintf(wr, "%24.16e %24.16e %24.16e\n", t, k*fstep, v);
          }
        fprintf(wr, "\n");
      }
    else
      { fprintf(stderr, "unspecified format, aborted\n");
        exit(1);
      }
  }

void write_output_trailer(FILE *wr, Output_Fmt_t ofmt)
  { 
    if (ofmt == OFM_PGM)
      { /* Nothing to do. */
      }
    else if (ofmt == OFM_ASCII)
      { /* Nothing to do. */
      }
    else
      { fprintf(stderr, "unspecified format, aborted\n");
        exit(1);
      }
    fflush(wr);
  }
  
void compute_frame
  ( jsaudio_t *si, 
    uint32_t skip, 
    Window_Type_t wtype,
    uint32_t nw,
    fftw_complex *in, 
    fftw_complex *out, 
    fftw_plan *plan,
    double pwr[]
  )
  {
    uint32_t nfrq = nw/2 + 1; /* Number of distinct frequencies. */

    /* Power spectrum normalization factor: */
    double norm = 1.0/nw;

    /* Extract samples {0..nw-1} of frame: */
    assert(skip + nw <= si->ns); 
    for (uint32_t k = 0;  k < nw; k++) { in[k][0] = si->sv[0][skip+k]; in[k][1] = 0.0; }
    
    /* Apply window function: */
    if (wtype == WIN_RECT) 
      { /* Nothing to do. */ }
    else if (wtype == WIN_HANN)
      { /* Multiply {tr} by the Hann window: */
        for (uint32_t k = 0;  k < nw; k++) 
          { double x = M_PI*(2*((double)k)/((double)nw) - 1.0);
            double w = 0.5*(1.0 + cos(x)); 
            in[k][0] *= w; in[k][1] *= w;
          }
      }
    
    /* Compute discrete transform: */
    fftw_execute(*plan);
    
    /* Save power spectrum in {pwr}: */
    for (uint32_t f = 0;  f < nfrq; f++)
      { double rep = out[f][0], imp = out[f][1];
        double sum = rep*rep + imp*imp;
        if ((f > 0) && (f < nw - f))
          { /* Term with frequency {-f}: */
            uint32_t g = nw - f;
            double rem = out[g][0], imm = out[g][1];
            sum += rem*rem + imm*imm;
          }
        pwr[f] = sum * norm;          
      } 
  }

void update_egram_range(uint32_t nfrq, double pwr[], double *emin, double *emax)
  {
    uint32_t f;
    for (f = 0; f < nfrq; f++)
      { double v = pwr[f]; 
        if (v < (*emin)) { (*emin) = v; }
        if (v > (*emax)) { (*emax) = v; }
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
    
    argparser_get_keyword(pp, "-input");
    if (argparser_keyword_present_next(pp, "ascii"))
      { o->ifmt = IFM_ASCII; }
    else if (argparser_keyword_present_next(pp, "au"))
      { o->ifmt = IFM_AU; }
    else
      { argparser_error(pp, "missing or invalid input format"); }
      
    argparser_get_keyword(pp, "-window");
    if (argparser_keyword_present_next(pp, "hann"))
      { o->wtype = WIN_HANN; }
    else if (argparser_keyword_present_next(pp, "rect"))
      { o->wtype = WIN_RECT; }
    else
      { argparser_error(pp, "missing or invalid window function"); }
    o->wsize = (uint32_t)argparser_get_next_int(pp, 1, MAX_WSIZE);
    
    if (argparser_keyword_present(pp, "-stride"))
      { o->stride = (uint32_t)argparser_get_next_int(pp, 0, MAXINT); }
    else 
      { o->stride = o->wsize/8; }

    argparser_get_keyword(pp, "-output");
    if (argparser_keyword_present_next(pp, "pgm"))
      { o->ofmt = OFM_PGM; }
    else if (argparser_keyword_present_next(pp, "ascii"))
      { o->ofmt = OFM_ASCII; }
    else
      { argparser_error(pp, "missing or invalid output format"); }
    
    if (argparser_keyword_present(pp, "-range"))
      { o->vmin = argparser_get_next_double(pp, 0.0, DBL_MAX);
        o->vmax = argparser_get_next_double(pp, o->vmin, DBL_MAX);
      }
    else 
      { o->vmin = 0.0;
        o->vmax = 1.0;
      }

    argparser_finish(pp);
    return o;
  }
