#define PROG_NAME "auclean"
#define PROG_DESC "cleanup noise from audio files"
#define PROG_VERS "1.0"

/* Copyright © 2006 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2023-03-02 19:31:31 by stolfi */

#define PROG_HELP \
  PROG_NAME "\\\n" \
  "    -window {WINDOW_SIZE}  [ -overlap {WINDOW_OVERLAP} ] \\\n" \
  "    [ -noiseFile {NOISE_FILE} .. ] \\\n" \
  "    [ -voice {VOICE_GAIN} \\\n" \
  "    [ -writeChaff ] \\\n" \
  "    [ -debugTime {DEBUG_TIME} ] \\\n" \
  "    -outPrefix {OUT_PREFIX} \\\n" \
  "    " argparser_help_info_HELP "\\\n" \
  "    < {INPUT_FILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  "  " PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads an audio file and outputs a cleaned-up version of it.\n" \
  "\n" \
  "  The input signal is read from standard input and is converted from" \
  " its external format to a sequence of double-precision floating-point" \
  " samples.  This sequence is processed in overlapping batches (frames)" \
  " of {WINDOW_SIZE} consecutive samples, where successive batches are" \
  " displaced by {WINDOW_SIZE}/{WINDOW_OVERLAP} samples.  The program" \
  " multiplies each" \
  " frame by a smooth partition-of-unity windowing function, computes" \
  " its discrete Fourier transform (DFT), discards unwanted components," \
  " reduces the modulus of each Fourier component to remove random" \
  " noise, computes the inverse transform, and" \
  " adds it to the output signal, which is written to the output file.\n" \
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
  "  -noiseFile {NOISE_FILE}\n" \
  "    Specifies the name of a \".au\" file containing a typical" \
  " example of the noise to be removed.  This parameter is optional" \
  " and may be repeated several times to provide examples of different" \
  " types of noise.\n" \
  "\n" \
  "  -voice {VOICE_GAIN} \n" \
  "    Tries to enhance voice-like sounds, recognized by a" \
  " fundamental frequency in the range 80 to 11100 Hz and" \
  " strong harmonics thereof.  The default is no" \
  " voice-like enhancement.\n" \
  "\n" \
  "  -writeChaff\n" \
  "    If present, this otional parameter specifies that the output must be the noise that would" \
  " be removed, instead of the cleaned-up signal. \n" \
  "\n" \
  "  -debugTime {DEBUG_TIME}\n" \
  "    If present, this otional parameter asks for debugging" \
  " output to be produced for the frame closest to the specified" \
  " time coordinate (in seconds). \n" \
  "\n" \
  argparser_help_info_HELP_INFO "\n" \
  "SEE ALSO\n" \
  "  spegram(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created oct/2006 by Jorge Stolfi, State University" \
  " of Campinas (UNICAMP)." \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-03-01 J.Stolfi Ouput prefix option instead of {stdout}." \
  "  2023-02-24 J.Stolfi Split off the lineat time-invariant filters to {aufilter.c}." \
  "  2023-02-24 J.Stolfi Multiple noise examples and robust LSQ fitting.\n" \
  "  2023-02-22 J.Stolfi Improved noise matching algorithm.\n" \

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
#include <jsfile.h>
#include <rn.h>
#include <rmxn.h>
#include <gauss_elim.h>
#include <argparser.h>
#include <jsrandom.h>

#include <jsaudio.h>
#include <jsaudio_au.h>

/* DATA TYPES */

/* Command line options: */
typedef struct options_t
  { /* Fourier analysis parameters: */
    int32_t wsize;          /* Window width. */
    int32_t overlap;        /* Num of consecutive overlapping frames. */
    /* Parameters of the noise filter: */
    string_vec_t noiseFile;  /* Names of files with noise examples (.au format) or NULL. */
    double voice;            /* Enhancement factor for voice-like sounds. */
    /* Output parameters: */
    char *outPrefix;         /* Prefix for output file names. */
    bool_t writeChaff;            /* TRUE to write the removed noise signal. */
    /* Debugging parameters: */
    double debugTime;     /* Time coordinate of frame to debug, or negative if none. */
  } options_t;
  
/* PROTOTYPES */

int32_t main(int32_t argc, char **argv);
  /* Main program. */

sound_t read_signal(char *fname, char *sname);
  /* Reads a sound clip from the file named {fname} (or from {stdin},
    if {fname} is "stdin"). The file must be in the Sun ".au" audio file format.
    The {sname} is the nature of the signal, e.g. "input" or "noise example". */

void write_signal(char *prefix, char *tag, sound_t *s);
  /* Writes the sound clip {s} to file "{prefix}-{tag}.au" in the Sun ".au" audio file format. */

int32_t get_noise_spectra
  ( int32_t nw, 
    int32_t no, 
    int32_t nnoi, 
    char *fname[], 
    double pwr[],
    /* FFTW work areas and plan: */
    fftw_complex in[],    /* FFTW input vector. */
    fftw_complex out[],   /* FFTW output vector. */
    fftw_plan *plan_dir   /* Precomputed direct FFTW parameters ({in ==> out}). */
  );
  /* Reads {nnoi} examples of noise from files whose names are {fname[0..nnoi-1]}, 
    computes their average windowed power spectra 
    and stores them into the rows of the matrix {pwr[0..nnoi-1,0..nfrq-1]}.
    Each noise spectrum is normalized to unit total power.
    
    Noise spectra with not sufficiently independent of the previous ones
    may be discarded. Returns the number of spectra (rows of {noise_pwr})
    actually retained. */

void process_frame
  ( sound_t *si,              /* Input sound signal. */
    int32_t skip,             /* Index of first sample of frame to be processed. */
    int32_t nw,               /* Samples per frame (window). */
    int32_t no,               /* Number of successive overlapping frames. */
    int32_t nnoi,             /* Number of noise examples. */
    double noise_pwr[],       /* Power spectra of noise examples. */
    double voice,             /* Enhancement factor for voice-like sounds. */
    sound_t *so,              /* Output: denoised sound signal. */
    sound_t *sn,              /* Output: removed noise. */
    bool_t debug,             /* Generate diagnostic output. */
    /* FFTW work areas and plan: */
    fftw_complex in[],        /* WORK: FFT input vector. */
    fftw_complex out[],       /* WORK: FFT output vector. */
    fftw_plan *plan_dir,      /* Precomputed direct FFTW parameters ({in ==> out}). */
    fftw_plan *plan_inv       /* Precomputed inverse FFTW parameters ({out ==> in}). */
  );
  /* Processes one frame of the input signal {si}, consisting of {nw}
    samples starting at sample {skip}, and pastes the filtered result
    onto the output signal {so}. Assumes that the frame overlaps the
    next {no-1} frames. 
    
    If {sn.ns} is not zero and large enough, also pastes
    onto {sn} the noise signal that was removed.
    
    The processing consists of applying a fixed filter with the
    transfer function {fixgain}, followed by a noise-removal 
    step, and a voice enhancement step. 
    
    The noise removal uses {nnoi} noise examples whose power spectra are
    assumed to be in the first {nnoi} rows of the matrix
    {noise_pwr[0..nnoi-1,0..nfrq-1]}.
    
    The FFT work vectors {in} and {out} must be allocated by the
    caller. */
    
void apply_hann_window_to_signal(double frame[], int32_t nw, int32_t no);
  /* Multiplies the signal {frame[0..nw-1]} by a Hann smoothing window,
    scaled so that the windows of all frames, assuming {no} overlap factor,
    add up to the unit constant function. */
  
void apply_filter
  ( int32_t nw,
    fftw_complex frame_ft[],
    double gain[],
    fftw_plan *plan_inv,
    fftw_complex in[],
    fftw_complex out[],
    double frout[]
  );
  /* Applies the filter {gain[0..nw-1]} to the frame's FT {frame_ft[0..nw-1]}
    and computes the inverse Fourier transform.  The result, which must be 
    real and apodized, is returned in {frout[0..nw-1]}. */
  
void splat_frame
  ( int32_t nw,
    double frame[],
    int32_t ic,
    int32_t skip,
    sound_t *sd
  );
  /* Adds the windowful of samples {frame[0..nw-1]} to channel {ic} of the 
    sound signal {sd}, shifted by {skip}. */

void estimate_noise_in_frame
  ( int32_t nw, 
    double frame_pwr[], 
    int32_t nnoi,
    double noise_pwr[],
    double audio_pwr[],
    double chaff_pwr[],
    double fsmp,
    bool_t debug
  );
  /* Tries to separate the power spectrum {frame_pwr[0..nfrq-1]} of a
    frame into the spectrum {audio_pwr[0..nfrq-1]} of a pure signal
    ("audio") plus a "chaff" spectrum {chaff_pwr[0..nfrq-1]}, the latter
    being a linear combination of {nnoi} "noise" signals whose average
    power spectra are the rows of the array
    {noise_pwr[0..nnoi-1,0..nfrq-1]}.
  
    The separation relies on the assumption that the 
    audio signal has practically zero power at some frequencies.
    If that is not the case, some of the audio signal
    may be misidentified as chaff.  */
    
    
void compute_freq_weights(int32_t nw, double frame_pwr[], double chaff_pwr[], double wt[]);
  /* Computes weights {wt[0..nfrq-1]} for each frequency for least squares 
    fitting of noise spectrum combination to frame spectrum.
    
    The weight {wt[kf]} is based on the frame's power spectrum {frame_pwr[kf]} and
    some estimate {chaff_pwr[kf]} of the noise power in that fequency. Only the 
    relative amounts of {chaff_pwr[0..nfrq-1]} matter.  The weight is higher
    where the chaff power is higher and the frame power is lower. 
    The weights are normalized to add to 1. */

void compute_power_spectrum(int32_t nw, fftw_complex S[], double pwr[]);
  /* Computes the power spectrum {pwr[0..nfrq-1]} from the Fourier transform
    {S[0..nw-1]} where {nfrq = nw/2 + 1}. */

void get_noise_spectrum
  ( char *fname, 
    int32_t nw, 
    int32_t no,
    double pwr[],
    /* FFTW work areas and plan: */
    fftw_complex in[],  /* FFTW input vector. */
    fftw_complex out[], /* FFTW output vector. */
    fftw_plan *plan_dir
  );
  /* Reads a sample signal from the ".au" file {fname} and computes its
    average power spectrum {pwr[0..nfrq-1]}, by using windows of with
    {nw} with {no}-fold overlap. The spectrum is normalized to have unit total power.
    The number of frequencies {nfrq} is, by
    Nyquist, {nw/2 + 1}. */
    
void show_filter(int32_t nw, double gain[]);
  /* Writes to {stderr} the transfer function {gain[0..nw]} of some filter. */

options_t* get_options(int32_t argc, char **argv);
  /* Parses the command-line arguments. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    /* Read the input signal: */
    sound_t si = read_signal("stdin", "input");
    
    int32_t nw = o->wsize;    /* Number of samples in window. */
    fprintf(stderr, "window size = %d samples (%.6f sec)\n", nw, ((double)nw)/si.fsmp);
    int32_t nfrq = nw/2 + 1;   /* Number of distinct frequencies, per Nyquist. */
    fprintf(stderr, "spectrum has %d frequencies {0..%d}\n", nfrq, nfrq-1);
    int32_t no = o->overlap;  /* Number of consecutive samples that overlap. */
    fprintf(stderr, "frame overlap count = %d\n", no);
    assert(no % 2 == 0);
    assert(nw % no == 0);
    int32_t stride = nw/no; /* Number of time steps between spectrograms. */
    int32_t nfrm = (si.ns - nw)/stride + 1;  /* Number of frames in input signal. */
    fprintf(stderr, "number of frames = %d\n", nfrm);
    fprintf(stderr, "frame spacing = %d samples (%.6f sec)\n", stride, ((double)stride)/si.fsmp);

    /* Allocate working storage for Fast Fourier Transform (FFT): */
    fftw_complex in[nw];  /* FFTW input vector. */
    fftw_complex out[nw]; /* FFTW output vector. */

    /* Precompute parameters for the FFT: */
    fftw_plan plan_dir = fftw_plan_dft_1d(nw, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_plan plan_inv = fftw_plan_dft_1d(nw, out, in, FFTW_BACKWARD, FFTW_ESTIMATE);

    /* Get the noise spectra {noise_pwr[0..nnoi-1,0..nfrq-1]}: */
    int32_t nnoi = o->noiseFile.ne; /* Number of noise clips provided. */
    double *noise_pwr = rmxn_alloc(nnoi, nfrq); /* Rows are noise spectra in log scale. */
    nnoi = get_noise_spectra(nw, no, nnoi, o->noiseFile.e, noise_pwr, in, out, &plan_dir);
    
    /* Allocate and and initialize the output signals: */
    sound_t so = jsa_allocate_sound(si.nc, si.ns);
    so.fsmp = si.fsmp;
    sound_t sn = jsa_allocate_sound(si.nc, (o->writeChaff ? si.ns : 0));
    sn.fsmp = si.fsmp;
    
    /* Compute frames and output them: */
    double totTime = si.ns/((double)si.fsmp);
    int32_t debug_ifrm = (int32_t)floor(nfrm*(o->debugTime/totTime));
    for (int32_t ifrm = 0; ifrm < nfrm; ifrm++)
      { /* Compute index {skip} of first sample in frame: */
        int32_t skip = ifrm*stride;
        bool_t debug = (ifrm == debug_ifrm);
        process_frame
          ( &si, skip, 
            nw, no, 
            nnoi, noise_pwr, 
            o->voice,
            &so,
            &sn,
            debug,
            in, out, &plan_dir, &plan_inv
          );
      }
      
    /* Output sound: */
    write_signal(o->outPrefix, "cln", &so);
    if (o->writeChaff) { write_signal(o->outPrefix, "nse", &sn); }

    return 0;
  }
  
sound_t read_signal(char *fname, char *sname)
  { fprintf(stderr, "Reading %s file from %s...\n", sname, fname);
    FILE *rd = (strcmp(fname, "stdin") == 0 ? stdin : open_read(fname, TRUE));
    sound_t si = jsa_read_au_file(rd);
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

int32_t get_noise_spectra
  ( int32_t nw, 
    int32_t no, 
    int32_t nnoi, 
    char *fname[], 
    double pwr[],
    /* FFTW work areas and plan: */
    fftw_complex in[],  /* FFTW input vector. */
    fftw_complex out[], /* FFTW output vector. */
    fftw_plan *plan_dir
  )
  { 
    int32_t nfrq = nw/2 + 1; 
    double npwrk[nfrq];
    int32_t m = 0; /* Number of noise examples actually retained. */
    for (int32_t i = 0; i < nnoi; i++)
      { /* Read the noise example {i} and compute the average spectrum of its frames: */
        fprintf(stderr, "reading noise example %d\n...", i);
        get_noise_spectrum(fname[i], nw, no, npwrk, in, out, plan_dir);
        for (int32_t kf = 0; kf < nfrq; kf++)
          { pwr[m*nfrq + kf] = npwrk[kf]; }
        m++;
      }
    fprintf(stderr, "retained %d of %d noise example spectra\n", m, nnoi); 
    return m;
  }

void get_noise_spectrum
  ( char *fname, 
    int32_t nw, 
    int32_t no,
    double pwr[],
    /* FFTW work areas and plan: */
    fftw_complex in[],  /* FFTW input vector. */
    fftw_complex out[], /* FFTW output vector. */
    fftw_plan *plan_dir
  )
  { sound_t sn = read_signal(fname, "noise example");
    /* Get displacement {stride} between frames: */
    int32_t stride = nw/no;
    /* Get number {nfrm} of frames in noise spectrum: */
    int32_t nfrm = (sn.ns - nw)/stride + 1; 
    fprintf(stderr, "number of frames in noise example = %d\n", nfrm);
    /* Get the number of channels: */
    int32_t nc = sn.nc;

    /* Get the number of distinct frequencies: */
    int32_t nfrq = nw/2 + 1;

    /* Initialize the total spectrum: */
    for (int32_t kf = 0; kf < nfrq; kf++) { pwr[kf] = 0; }

    /* Compute frame spectrograms of noise signal and accumulate them: */
    double pwri[nfrq]; /* Power spectrum of one frame. */
    for (int32_t ifrm = 0; ifrm < nfrm; ifrm++)
      { /* Compute index {skip} of first sample in frame: */
        int32_t skip = ifrm*stride;
     
        /* Loop on channels: */
        for (int32_t ic = 0; ic < nc; ic++)
          { 
            /* Extract frame from input signal: */
            assert(skip + nw <= sn.ns); 
            double frame[nw];           /* Windowful of samples from an input channel. */
            for (int32_t kw = 0; kw < nw; kw++) 
              { frame[kw] = sn.sv[ic][skip+kw]; }

            /* Apply the Hann unit-partition window function: */
            apply_hann_window_to_signal(frame, nw, no); 

            /* Perform Fourier transform: */
            for (int32_t kw = 0; kw < nw; kw++) 
              { in[kw][0] = frame[kw]; in[kw][1] = 0.0; }
            fftw_execute(*plan_dir);
            fftw_complex frame_ft[nw];  /* FR of {frame}. */
            for (int32_t kw = 0; kw < nw; kw++)
              { frame_ft[kw][0] = out[kw][0]; frame_ft[kw][1] = out[kw][1]; }
            
            /* Convert to power spectrum and accumulate: */
            compute_power_spectrum(nw, frame_ft, pwri);
            for (int32_t kf = 0; kf < nfrq; kf++)
              { pwr[kf] += pwri[kf]; }
          }
      }
      
    /* Convert sum to average: */
    double tot_pwr = rn_sum(nfrq, pwr) + 1.0e-200;
    for (int32_t kf = 0; kf < nfrq; kf++) 
      { assert(isfinite(pwr[kf]));
        pwr[kf] /= tot_pwr; 
        assert(isfinite(pwr[kf]));
      }
  }

void write_signal(char *prefix, char *tag, sound_t *s)
  { char *fname = NULL;
    asprintf(&fname, "%s-%s.au", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    fprintf(stderr, "Writing output file");
    fprintf(stderr, " (Sun \".au\" format)\n");
    jsa_write_au_file(wr, s);
    free(fname);
  }
  
void process_frame
  ( sound_t *si,              /* Input sound signal. */
    int32_t skip,             /* Index of first sample of frame to be processed. */
    int32_t nw,               /* Samples per frame (window). */
    int32_t no,               /* Number of successive overlapping frames. */
    int32_t nnoi,             /* Number of noise examples. */
    double noise_pwr[],       /* Power spectra of noise examples. */
    double voice,             /* Enhancement factor for voice-like sounds. */
    sound_t *so,              /* Output: denoised sound signal. */
    sound_t *sn,              /* Output: removed noise signal. */
    bool_t debug,             /* Generate diagnostic output. */
    /* FFTW work areas and plans: */
    fftw_complex in[],        /* WORK: FFT input vector. */
    fftw_complex out[],       /* WORK: FFT output vector. */
    fftw_plan *plan_dir,      /* Precomputed direct FFTW parameters ({in ==> out}). */
    fftw_plan *plan_inv       /* Precomputed inverse FFTW parameters ({out ==> in}). */
  )
  { 
    assert(si->fsmp == so->fsmp);
    int32_t nfrq = nw/2 + 1;
   
    /* Number of channels: */
    int32_t nc = (si->nc < so->nc ? si->nc : so->nc);
    
    /* Loop on channels: */
    for (int32_t ic = 0; ic < nc; ic++)
      { /* Extract frame from input signal: */
        assert(skip + nw <= si->ns); 
        double frame[nw];           /* Windowful of samples from an input channel. */
        for (int32_t kw = 0; kw < nw; kw++) 
          { frame[kw] = si->sv[ic][skip+kw]; }
        
        /* Apply the Hann unit-partition window function: */
        apply_hann_window_to_signal(frame, nw, no); 
      
        /* Perform Fourier transform: */
        for (int32_t kw = 0; kw < nw; kw++) 
          { in[kw][0] = frame[kw]; in[kw][1] = 0.0; }
        fftw_execute(*plan_dir);
        fftw_complex frame_ft[nw];  /* FR of {frame}. */
        for (int32_t kw = 0; kw < nw; kw++)
          { frame_ft[kw][0] = out[kw][0]; frame_ft[kw][1] = out[kw][1]; }
        
        /* Compute frame power spectrum: */
        double frame_pwr[nfrq];     /* Power spectrum of {frame}. */
        compute_power_spectrum(nw, frame_ft, frame_pwr);

        /* Partition {frame_pwr} into {audio_pwr} and {chaff_pwr}: */
        double audio_pwr[nfrq]; /* Noise removed. */
        double chaff_pwr[nfrq]; /* Noise removed. */
        estimate_noise_in_frame(nw, frame_pwr, nnoi, noise_pwr, audio_pwr, chaff_pwr, si->fsmp, debug);

        /* Compute denoising filter gain: */
        double gain[nw];
        for (int32_t kf = 0; kf < nfrq; kf++) 
          { int32_t jf = (nw - kf) % nw; /* The freq {-kf} mod {nw}. */
            double gk = sqrt(audio_pwr[kf]/frame_pwr[kf]);
            gain[kf] = gk;
            if (jf != kf) { gain[jf] = gk; }
          }
        
        /* Apply the filter to the frame's FT and invert: */
        double audio[nw];
        apply_filter(nw, frame_ft, gain, plan_inv, in, out, audio);
        splat_frame(nw, audio, ic, skip, so);
        
        if (sn->ns != 0)
          { double chaff[nw];
            for (int32_t kw = 0; kw < nw; kw++)
              { chaff[kw] = frame[kw] - audio[kw]; }
            splat_frame(nw, chaff, ic, skip, sn);
          }
      }
  }
          
void apply_filter
  ( int32_t nw,
    fftw_complex frame_ft[],
    double gain[],
    fftw_plan *plan_inv,
    fftw_complex in[],
    fftw_complex out[],
    double frout[]
  )
  {
    /* Apply the filter to the input signal, pepare for FT inv: */
    for (int32_t kw = 0; kw < nw; kw++) 
      { double gn = gain[kw];
        out[kw][0] = gn*frame_ft[kw][0]; 
        out[kw][1] = gn*frame_ft[kw][1];
      }

    /* Compute inverse FFT: */
    fftw_execute(*plan_inv);
    for (int32_t kw = 0; kw < nw; kw++) 
      { if(fabs(in[kw][1]) > 1.0e-6)
          { fprintf(stderr, "%4d %16.10f %16.10f\n", kw, in[kw][0], in[kw][1]);
            demand(FALSE, "filtered signal is complex");
          } 
        frout[kw] = in[kw][0]/nw;
      }
  }

void splat_frame
  ( int32_t nw,
    double frame[],
    int32_t ic,
    int32_t skip,
    sound_t *sd
  )
  { assert(skip + nw <= sd->ns); 
    for (int32_t kw = 0; kw < nw; kw++) 
      { sd->sv[ic][skip+kw] += frame[kw]; }
  }

void show_filter(int32_t nw, double gain[])
  { for (int32_t kw = 0; kw < nw; kw++) 
      { fprintf(stderr, "gain[%03d] = %16.12f\n", kw, gain[kw]); }
  }

void apply_hann_window_to_signal(double frame[], int32_t nw, int32_t no)
  {
    /* Scaling factor for Hann unit-partition window function */
    assert(nw % no == 0);
    double fac = 1.0/((double)no);

    for (int32_t kw = 0; kw < nw; kw++) 
      { double x = M_PI*(2*((double)kw)/((double)nw) - 1.0);
        double w = fac*(1.0 + cos(x)); 
        frame[kw] *= w; 
      }
  }
void compute_power_spectrum(int32_t nw, fftw_complex S[], double pwr[])
  {
    /* Number of frequencies: */
    int32_t nfrq = nw/2 + 1;
    
    /* Normalization factor for power spectrum: */
    double norm = 1.0/nw;
  
    for (int32_t kf = 0; kf < nfrq; kf++)
      { /* Get the freq {jf} equivalent to freq {-kf} mod {nw}. */
        int32_t jf = (nw - kf) % nw; 
        
        /* Get real and imaginay parts of coeffs of freqs {kf} and {-kf}: */
        double rek = S[kf][0], imk = S[kf][1];
        double pk = rek*rek + imk*imk;
        if (jf != kf) 
          { double rej = S[jf][0], imj = S[jf][1];
            pk += rej*rej + imj*imj;
          }
        assert(isfinite(pk));
        assert(pk >= 0);
        pwr[kf] = norm*pk;
      }
  }

void estimate_noise_in_frame
  ( int32_t nw, 
    double frame_pwr[], 
    int32_t nnoi,
    double noise_pwr[],
    double audio_pwr[],
    double chaff_pwr[],
    double fsmp,
    bool_t debug
  )
  { 
    if (debug) { fprintf(stderr, "starting {compose_denoising_filter} nnoi = %d\n", nnoi); }

    int32_t nfrq = nw/2 + 1; /* Number of distinct frequencies: */

    if (nnoi == 0)
      { /* Nothing to do: */
        for (int32_t kf = 0; kf < nfrq; kf++)
          { audio_pwr[kf] = frame_pwr[kf];
            chaff_pwr[kf] = 0;
          }
        return;
      }

    for (int32_t kf = 0; kf < nfrq; kf++)
      { /* Initialize {chaff_pwr} with an arbitrary multiple of est noise power for freq {kf}: */
        double sum_p2 = 1.0e-200;
        for (int32_t i = 0; i < nnoi; i++)
          { double pk = noise_pwr[i*nfrq + kf];
            sum_p2 += pk*pk;
          }
        chaff_pwr[kf] = sqrt(sum_p2);
        assert(chaff_pwr[kf] > 0);
      }

    /* The audio and chaff are identified by least squares fitting,
      with weights {wt[0..nfrq-1]} that emphasize frequencies
      where the audio power seems closer to noise power.  These
      weights are adjusted iteratively later. */
    double wt[nfrq];
    compute_freq_weights(nw, frame_pwr, chaff_pwr, wt);
  
    int32_t nnoi_u = nnoi; /* Number of noise samples to use in fit. */
    int32_t ixnoi[nnoi]; /* Indices of noise samples to use in fit. */
    for (int32_t i = 0; i < nnoi; i++) { ixnoi[i] = i; }

    int32_t maxiter = nnoi+1;
    int32_t niter = 0;
    while ((niter < maxiter) && (nnoi_u > 0))
      { niter++;
        if (debug) { fprintf(stderr, "iteration %d\n", niter); } 

        FILE *wr = NULL;
        if (debug)
          { char *fname = NULL;
            asprintf(&fname, "out/debug_%02d.pwr", niter);
            wr = open_write(fname, TRUE);
            free(fname);
          }

        /* Fit a combination of noise spectra to {frame_pwr} with weights {wt}. */
        /* The noise spectra do use are {ixnoi[0..nnoi_u-1]}. */

        /* Compute the moment matrix {M} and the indep vector {b}: */
        double M[nnoi_u*nnoi_u];
        double b[nnoi_u];
        for (int32_t i_u = 0; i_u < nnoi_u; i_u++)
          { int32_t i = ixnoi[i_u];
            double *npi = &(noise_pwr[i*nfrq]);
            double bi = 0;
            for (int32_t kf = 0; kf < nfrq; kf++)
              { bi += wt[kf]*npi[kf]*frame_pwr[kf]; }
            b[i_u] = bi;
            
            double *Mi = &(M[i_u*nnoi_u]);
            for (int32_t j_u = 0; j_u < nnoi_u; j_u++)
              { int32_t j = ixnoi[j_u];
                double *npj = &(noise_pwr[j*nfrq]);
                double Mij = 0;
                for (int32_t kf = 0; kf < nfrq; kf++)
                  { Mij += wt[kf]*npi[kf]*npj[kf]; }
                Mi[j_u] = Mij;
                assert(isfinite(Mij));
              }
          }
        /* Solve the system {M*x = b}: */
        double x[nnoi_u];
        int32_t rank = gsel_solve(nnoi_u, nnoi_u, M, 1, b, x, 1.0e-10);
        if ((debug) || (! isfinite(x[0])))
          { gsel_print_system(stderr, "%12.8f", "linear system:", nnoi_u, nnoi_u, M, 1, b, NULL);
            rn_gen_print(stderr, nnoi_u, x, "%12.8f", "solution: ", " ", "\n"); 
          }
        assert(isfinite(x[0]));
        if (rank < nnoi_u)
          { fprintf(stderr, "lin sys solution failed, rank = %d\n", rank); }
          
        /* Exclude any noise component that got zero or negative {x}: */
        int32_t nnoi_new = 0;
        for (int32_t i_u = 0; i_u < nnoi_u; i_u++)
          { if (! isfinite(x[i_u]))
              { fprintf(stderr, "!? x[%d] = %14.8f\n", i_u, x[i_u]);
                assert(FALSE);
              }
            if (x[i_u] > 0.0)
              { ixnoi[nnoi_new] = ixnoi[i_u]; 
                x[nnoi_new] = x[i_u];
                nnoi_new++;
              }
            else if (debug)
              { fprintf(stderr, "noise source %d", ixnoi[i_u]); 
                fprintf(stderr, "  coeff x[%d] = %14.10f, excluded\n", i_u, x[i_u]); 
              }
          }
        if (nnoi_new < nnoi_u)
          { maxiter += nnoi_u - nnoi_new; /* Do extra work if we excluded some noise samples. */
            nnoi_u = nnoi_new; 
          }
          
        /* Compute chaff as determined by solution of lin sys: */
        for (int32_t kf = 0; kf < nfrq; kf++)
          { double chk = 0;
            for (int32_t i_u = 0; i_u < nnoi_u; i_u++)
              { int32_t i = ixnoi[i_u];
                chk += x[i_u]*noise_pwr[i*nfrq + kf];
              }
            chk = fmax(1.0e-200, chk);
              
            /*  Compute the "audio" power: */
            double auk = fmax(1.0e-200, frame_pwr[kf] - 1.1*chk);
    
            if (debug)
              { double freq = (kf*fsmp)/nw;
                fprintf(wr, "%10.4f", freq);
                fprintf(wr, "  %14.8f %14.8f", frame_pwr[kf], chaff_pwr[kf]);
                fprintf(wr, "  %14.10f %14.8f %14.8f\n", wt[kf], chk, auk);
              }          

            assert(isfinite(chk));
            assert(isfinite(auk));
          
            chaff_pwr[kf] = chk;
            audio_pwr[kf] = auk;
          }

         if (debug) { fclose(wr); }      

        /* Recompute weights using {chaff_pwr} instead of {mush_pwr}: */
        compute_freq_weights(nw, frame_pwr, chaff_pwr, wt);
      }
  }
    
void compute_freq_weights(int32_t nw, double frame_pwr[], double chaff_pwr[], double wt[])
  { int32_t nfrq = nw/2 + 1;
    for (int32_t kf = 0; kf < nfrq; kf++)
      { double fpk = frame_pwr[kf];
        assert(fpk >= 0);
        double cpk = chaff_pwr[kf];
        assert(cpk >= 0);
        double rk = cpk*cpk/(fpk*fpk -2*fpk*cpk + 2*cpk*cpk + 1.0e-200);
        wt[kf] = cpk*rk;
        assert(isfinite(wt[kf]));
        assert(wt[kf] >= 0);
      }
    /* Normalize weights to unit sum: */
    double sum_wt = rn_sum(nfrq, wt) + 1.0e-200;
    rn_scale(nfrq, 1/sum_wt, wt, wt);
  }
            
#define MAX_WSIZE 9192

options_t* get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    argparser_get_keyword(pp, "-window");
    o->wsize = (int32_t)argparser_get_next_int(pp, 1, MAX_WSIZE);
    
    if (argparser_keyword_present(pp, "-overlap"))
      { o->overlap = (int32_t)argparser_get_next_int(pp, 2, MAXINT); }
    else 
      { o->overlap = 8; }
    if (o->overlap % 2 != 0) 
      { argparser_error(pp, "overlap count must be even"); }
    if (o->wsize % o->overlap != 0) 
      { argparser_error(pp, "overlap count must divide frame size"); }
    
    o->noiseFile = string_vec_new(10);
    int32_t nnoi = 0;
    while (argparser_keyword_present(pp, "-noiseFile"))
      { char *fn = argparser_get_next_non_keyword(pp);
        string_vec_expand(&(o->noiseFile), nnoi);
        o->noiseFile.e[nnoi] = fn;
        nnoi++;
      }
    string_vec_trim(&(o->noiseFile), nnoi);
    
    if (argparser_keyword_present(pp, "-voice"))
      { o->voice = argparser_get_next_double(pp, 0.0, DBL_MAX); }
    else 
      { o->voice = 0.0; }
    
    o->writeChaff = argparser_keyword_present(pp, "-writeChaff");
    
    if (argparser_keyword_present(pp, "-debugTime"))
      { o->debugTime = argparser_get_next_double(pp, -DBL_MAX, DBL_MAX); }
    else 
      { o->debugTime = -1.0; }
    
    argparser_get_keyword(pp, "-outPrefix");
    o->outPrefix = argparser_get_next_non_keyword(pp);
    
    argparser_finish(pp);
    return o;
  }
