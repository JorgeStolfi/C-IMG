#define PROG_NAME "pnmxhist"
#define PROG_DESC "compute a weighted histogram of pixel values in a PNM image"
#define PROG_VERS "1.0"

#define pnmxhist_C_COPYRIGHT \
  "Copyright © 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2025-08-02 12:14:04 by stolfi */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [ -channel {CHANNEL} ] \\\n" \
  "  [ -ignore {VALUE} ]... \\\n" \
  "  [ -mask {MASKFILE} ] \\\n" \
  "  [ -verbose ] \\\n" \
  "  [ INFILE ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads a PBM, PGM, or PPM file {INFILE} and writes to standard output" \
  " a histogram of its sample values.  If {INFILE} is omitted or \"-\", reads" \
  " the image from standard input.\n" \
  "\n" \
  "  In the case of a PPM" \
  " file, the histogram is built for a user-specified channel only.\n" \
  "\n" \
  "  The user may also specify a monochromatic weight (opacity) mask and/or a list of sample" \
  " values to be ignored when counting.  If the mask is given, its" \
  " samples are interpreted as fractional weights between 0 and 1; and each" \
  " sample of the input image whose weight is {W} is counted " \
  " as {W} pixels instead of one.\n" \
  "DESCRIPTION\n" \
  "  The output has seven columns\n" \
  "\n" \
  "    {VALUE} {COUNT} {COUNT_REL}  {ACCUM} {ACCUM_REL}  {MUCCA} {MUCCA_REL}\n" \
  "\n" \
  "  where\n" \
  "\n" \
  "    {VALUE} is a sample value ranging from 0 to the input image's {MAXVAL};\n" \
  "    {COUNT} is the number (or total weight) of samples with that value;\n" \
  "    {COUNT_REL} is {COUNT} divided by the total {COUNT} of all values;\n" \
  "    {ACCUM} is the total {COUNT}s of sample values less than {VALUE};\n" \
  "    {ACCUM_REL} is {ACCUM} divided by the total {COUNT};\n" \
  "    {MUCCA} is the total {COUNT}s of sample values greater than {VALUE};\n" \
  "    {MUCCA_REL} is {ACCUM} divided by the total {COUNT}.\n" \
  "\n" \
  "  Sample values which were specified through the \"-ignore\" option" \
  " are nor counted, but are listed in the output, always with zero {COUNT}.\n" \
  "\n" \
  "  The fields {COUNT}, {ACCUM}, and {MUCCA} may be fractional, and are" \
  " formatted with a sufficient number of decimal digits after the point.  The" \
  " fields {ACCUM} and {MUCCA} fields include only one half of the {COUNT} of the same line.\n" \
  "\n" \
  "  If the sum of all {COUNT}s is zero, the fractions {ACCUM_REL} and" \
  " {MUCCA_REL} are computed as if the counts of non-ignored values were equal" \
  " to some small value.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -channel {CHANNEL}\n" \
  "    This argument is alllowed and mandatory only for PPM input, and specifies" \
  " the channel to be examined.  The {CHANNEL} may be 0, \"R\", or \"red\";" \
  " 1, \"G\", \"grn\", or \"green\"; or 2, \"B\", \"blu\", or \"blue\".\n" \
  "\n" \
  "  -ignore {VALUE}\n" \
  "    This optional argument specifies that samples with the given {VALUE}" \
  " should be ignored when building the histogram.  The {VALUE} must" \
  " be in the range from 0 to the input image's {MAXVAL}.  This option" \
  " may be repeated in order to exclude two or more values.\n" \
  "\n" \
  "  -mask {MASKFILE}\n" \
  "    This optional argument specifies a PGM or PBM file that defines the weight of" \
  " each sample in the selected channel. Each of its samples is linearly mapped" \
  " from the range {0..MAXVAL} to a weight in the range {[0_1]}.  A pixel with" \
  " weight {W} contributes {W} to the respective {COUNT} field. If the mask is" \
  " not given, all non-ignored input samples are assumed to have weight 1.\n" \
  "\n" \
  "  -verbose\n" \
  "    This option causes the the program to print some informative messages.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgmhist(1), ppmhist(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created on 2010-08-08 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2023-03-07 by J. Stolfi: fixed channel selection bug. General cleanup.\n" \
  "  2010-08-08 by J. Stolfi: created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmxhist_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h> 
#include <jspnm.h> 
#include <jsfile.h> 
#include <jsmath.h> 
#include <uint16_image.h> 
#include <argparser.h> 

typedef struct options_t 
  { int32_t channel;         /* Channel to process, or -1 if not specified. */
    bool_vec_t ignore;   /* {ignore.e[v]} is true to ignore samples with value {v}. */
    bool_t verbose;      /* TRUE to mumble while working. */
    char *wt_filename;   /* Name of mask file, of NULL if not given. */
    char *in_filename;   /* Name of input file, or "-" if not given. */
  } options_t;
  
typedef struct output_format_t
  { int32_t sz_value;                    /* {VALUE}. */
    int32_t sz_count_abs, pr_count_abs;  /* {COUNT} */
    int32_t sz_count_rel, pr_count_rel;  /* {COUNT_REL} */
    int32_t sz_accum_abs, pr_accum_abs;  /* {ACCUM,MUCCA} */
    int32_t sz_accum_rel, pr_accum_rel;  /* {ACCUM_REL,MUCCA_REL} */
  } output_format_t;
  /* Specifies the width and precision of output fields.  */

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);
options_t *get_options(int32_t argc, char **argv);

void print_histogram
  ( FILE *wr,
    uint64_t count[],
    uint16_t in_maxval,
    uint16_t wt_maxval,
    bool_vec_t ignore,
    bool_t verbose
  );
  /* Prints the histogram {count[0..in_maxval]} to file {wr}. Assumes
    that that a sample value {v} is to be ignored if {ignore.e[v]} is
    true. Otherwise assumes that {count[v]} is the count of
    occurrences of value {v} in the input image, multiplied by the
    respective raw weights which range in {0..wt_maxval}.
    Therefore {count[v]} must be divided by {wt_maxval} to give the 
    actual weighted count. */ 

void choose_output_format
  ( uint16_t in_maxval,  /* {MAXVAL} of input image. */
    uint16_t wt_maxval,  /* {MAXVAL} of weight mask image, or 1 if none. */
    int32_t num_values,          /* Number of non-ignored input sample values. */
    uint64_t tot_wt_count,   /* Total weight of all non-ignored pixels. */
    uint64_t max_wt_count,   /* Maximum weight sum among all input sample values. */
    output_format_t *fmt     /* (OUT) The chosen output format parameters. */
  );
  /* Chooses the histogram output format parameters based on the
    given histogram parameters. */
  
/* ROUTINES */

#define MAX_WT_COUNT (((1LLU) << 50) - 1)
  /* The maximum allowed weighted count.  Must be exactly
    representable in a {double}, with one spare bit. */
 
int32_t main(int32_t argc, char **argv)
  {
    options_t *o = get_options(argc, argv);
    
    /* Open input image and read header: */
    FILE *in_file = open_read(o->in_filename, o->verbose);
    uint32_t chns, cols, rows;
    bool_t in_raw, in_bits;
    uint16_t in_maxval;
    { pnm_format_t in_format;
      pnm_read_header(in_file, &cols, &rows, &chns, &in_maxval, &in_raw, &in_bits, &in_format);
    }
    uint64_t pixs = (uint64_t)(cols*rows);
    if (o->verbose)
      { fprintf(stderr, "input image has %d colums, %d rows, %ld pixels\n", cols, rows, pixs); 
        fprintf(stderr, "input samples range in 0..%d\n", in_maxval);
      }
    
    /* Open weight image, if any: */
    FILE *wt_file = NULL;
    uint16_t wt_maxval;
    bool_t wt_raw, wt_bits;
    if (o->wt_filename != NULL)
      { /* Open weight image, and read header: */
        wt_file = open_read(o->wt_filename, o->verbose);
        uint32_t wt_chns, wt_cols, wt_rows;
        pnm_format_t wt_format;
        pnm_read_header(wt_file, &wt_cols, &wt_rows, &wt_chns, &wt_maxval, &wt_raw, &wt_bits, &wt_format);
        if (wt_chns != 1)
          { pnm_error("mask image must be monochromatic"); }
        if ((wt_cols != cols) || (wt_rows != rows))
          { pnm_error("mask image has incompatible size"); }
        if (o->verbose)
          { fprintf(stderr, "weight mask samples range in 0..%d\n", wt_maxval); }
      }
    else
      { wt_maxval = 1;        /* Pretend that the mask is a bitmap of all 1s: */
        wt_raw = wt_bits = 0; /* Just in case. */
      }
      
    /* Check for possible overflow: */
    if (((MAX_WT_COUNT / pixs) / (uint64_t)wt_maxval) == 0)
      { pnm_error("too many pixels and/or too large weights, counts may overflow"); }
    
    /* Check whether the "-ignore" values are all valid: */
    { for (int32_t v = in_maxval+1; v < o->ignore.ne; v++)
        { if (o->ignore.e[v]) 
            { pnm_error("value to ignore %llu is not valid", v); }
        }
    }
    
    /* Get the channel to operate on: */
    int32_t chn;
    if (chns == 1)
      { if (o->channel >= 0)
          { pnm_error("option \"-channel\" not allowed on monochromatic images"); }
        chn = 0;
      }
    else if (chns == 3)
      { if (o->channel < 0)
          { pnm_error("option \"-channel\" is required for color images"); chn = 0; }
        else
          { chn = o->channel; }
      }
      
    /* Allocate sample buffers for reading: */
    uint16_t *in_pix = uint16_image_alloc_pixel_row(cols, chns);
    uint16_t *wt_pix = uint16_image_alloc_pixel_row(cols, chns);

    /* Allocate histogram: */
    uint64_t count[in_maxval+1];
    for (int32_t v = 0;  v <= in_maxval; v++) { count[v] = 0; }
    
    /* Count pixels: */
    for (int32_t y = 0;  y < rows; y++)
      { pnm_read_pixels(in_file, in_pix, cols, chns, in_maxval, in_raw, in_bits);
        if (wt_file != NULL)
          { pnm_read_pixels(wt_file, wt_pix, cols, 1, wt_maxval, wt_raw, wt_bits); }
        for (int32_t x = 0;  x < cols; x++)
          { int32_t vxy = in_pix[x*((int32_t)chns) + chn];
            int32_t wxy = (wt_file == NULL ? 1 : wt_pix[x]);
            if (o->ignore.e[vxy]) { wxy = 0; }
            count[vxy] += (uint64_t)wxy;
          }
      } 
    
    /* Close input files: */
    if (wt_file != NULL) { fclose(wt_file); }
    if (in_file != stdin) { fclose(in_file); }
    
    print_histogram(stdout, count, in_maxval, wt_maxval, o->ignore, o->verbose);

    return(0);
  } 

void print_histogram
  ( FILE *wr,
    uint64_t count[],
    uint16_t in_maxval,
    uint16_t wt_maxval,
    bool_vec_t ignore,
    bool_t verbose
  )
  {
    /* Compute total and maximum weighted counts, and num of non-ignored values: */
    int32_t num_values = 0; /* Number of non-ignored sample values. */
    int32_t num_values_present = 0; /* Num of sample values with nonzero count. */
    int32_t min_value_present = in_maxval+1; /* Smallest sample value with nonzero count. */
    int32_t max_value_present = -1;          /* Largest sample value with nonzero count. */
    uint64_t tot_wt_count = 0; /* Sum of all {count}s. */
    uint64_t max_wt_count = 0; /* Max {count} entry. */
    uint16_t most_pop_value = 0; /* Sample value with hightest count. */
    for (int32_t v = 0;  v <= in_maxval; v++)
      { tot_wt_count += count[v];
        if (count[v] > max_wt_count) { max_wt_count = count[v]; most_pop_value = (uint16_t)v; }
        if (! ignore.e[v]) { num_values++; }
        if (count[v] != 0)
          { num_values_present++;
            if (min_value_present > in_maxval) { min_value_present = v; }
            max_value_present = v;
          }
      }
    if (num_values == 0)
      { pnm_error("cannot ignore *all* input values"); }
    if (verbose)
      { fprintf(stderr, "tot weighted count  %17.6f\n", ((double)tot_wt_count)/((double)wt_maxval));
        fprintf(stderr, "max weighted count  %17.6f", ((double)max_wt_count)/((double)wt_maxval));
        fprintf(stderr, " for sample value %d\n", most_pop_value);
        fprintf(stderr, "number valid sample values     %6d\n", num_values);
        if (tot_wt_count > 0)
          { 
            fprintf(stderr, "distinct valid values present %6d\n", num_values_present);
            fprintf(stderr, "actual sample range %d..%d\n", min_value_present, max_value_present);
          }
      }
    
    /* Determine the output field widths and precision: */
    output_format_t fmt;
    choose_output_format(in_maxval, wt_maxval, num_values, tot_wt_count, max_wt_count, &fmt);

    /* Write histogram: */
    assert(num_values > 0);
    int32_t lss_num = 0;           /* Number of smaller non-ignored values. */
    int32_t gtr_num = num_values;  /* Number of higher non-ignored values. */
    uint64_t lss_wt_count = 0;             /* Sum of {count}s for all smaller values */
    uint64_t gtr_wt_count = tot_wt_count;  /* Sum of {count}s for all higher values. */
    for (uint32_t v = 0;  v <= in_maxval; v++) 
      { if (ignore.e[v]) { assert(count[v] == 0); }
        /* Update {gtr_num,gtr_wt_count}: */
        if (! ignore.e[v]) { gtr_num--; }
        gtr_wt_count -= count[v];
        /* Compute the histogram fields: */
        double ct_abs = ((double)count[v])/((double)wt_maxval);
        double ac_abs = ((double)lss_wt_count + 0.5*(double)count[v])/((double)wt_maxval);
        double mu_abs = ((double)gtr_wt_count + 0.5*(double)count[v])/((double)wt_maxval);
        double ct_rel, ac_rel, mu_rel;
        if (tot_wt_count == 0)
          { ct_rel = 0;
            ac_rel = ((double)lss_num + 0.5)/((double) num_values);
            mu_rel = ((double)gtr_num + 0.5)/((double) num_values);
          }
        else
          { ct_rel = ((double)count[v])/((double)tot_wt_count);
            ac_rel = ((double)lss_wt_count + 0.5*(double)count[v])/((double)tot_wt_count);
            mu_rel = ((double)gtr_wt_count + 0.5*(double)count[v])/((double)tot_wt_count);
          }
        
        /* Print line: */
        fprintf(wr, "%*d", fmt.sz_value, (int32_t)v);
        fprintf(wr, " %*.*f", fmt.sz_count_abs, fmt.pr_count_abs, ct_abs);
        fprintf(wr, " %*.*f", fmt.sz_count_rel, fmt.pr_count_rel, ct_rel);
        fprintf(wr, " %*.*f", fmt.sz_accum_abs, fmt.pr_accum_abs, ac_abs);
        fprintf(wr, " %*.*f", fmt.sz_accum_rel, fmt.pr_accum_rel, ac_rel);
        fprintf(wr, " %*.*f", fmt.sz_accum_abs, fmt.pr_accum_abs, mu_abs);
        fprintf(wr, " %*.*f", fmt.sz_accum_rel, fmt.pr_accum_rel, mu_rel);
        fprintf(wr, "\n");

        /* Update {lss_num,lss_wt_count} for next iteration: */
        lss_wt_count += count[v];
        if (! ignore.e[v]) { lss_num++; }
      }
    fflush(wr);

    assert(lss_num == num_values);
    assert(gtr_num == 0);
    assert(lss_wt_count == tot_wt_count);
    assert(gtr_wt_count == 0);
  }
  
void choose_output_format
  ( uint16_t in_maxval,
    uint16_t wt_maxval,
    int32_t num_values, 
    uint64_t tot_wt_count, 
    uint64_t max_wt_count, 
    output_format_t *fmt
  )
  {
    /* Determine the output field width for {VALUE}: */
    fmt->sz_value = (int32_t)digits(in_maxval);

    /* Determine the output precision and field width for {COUNT}: */
    uint64_t max_count = (max_wt_count / wt_maxval); /* Max count, truncated to int32_t. */
    fmt->pr_count_abs = ((wt_maxval == 1) || (tot_wt_count == 0) ? 0 : (int32_t)digits(wt_maxval - 1)); 
    fmt->sz_count_abs = (int32_t)digits(max_count) + (int32_t)(fmt->pr_count_abs > 0) + (int32_t)fmt->pr_count_abs; 
    
    /* Determine the output precision and field width for {COUNT_REL}: */
    fmt->pr_count_rel = (tot_wt_count == max_wt_count ? 0 : (int32_t)digits(tot_wt_count - 1));
    fmt->sz_count_rel = 1 + (int32_t)(fmt->pr_count_rel > 0) + fmt->pr_count_rel;
    
    /* Determine the output precision and field width for {ACCUM,MUCCA}: */
    /* The cumulative counts are half-integers divided by {wt_maxval}: */
    uint64_t tot_count = (tot_wt_count / wt_maxval); /* Total count, truncated to int32_t. */
    fmt->pr_accum_abs = (wt_maxval == 1 ? 1 : (int32_t)digits((uint64_t)(2*wt_maxval - 1)));
    fmt->sz_accum_abs = (int32_t)digits(tot_count) + (int32_t)(fmt->pr_accum_abs > 0) + fmt->pr_accum_abs;
    
    /* Determine the output precision and field width for {ACCUM_REL,MUCCA_REL}: */
    if (tot_wt_count == 0)
      { /* The cumulative counts are half-integers divided by {num_values}: */
        fmt->pr_accum_rel = (num_values <= 1 ? 1 : (int32_t)digits((uint64_t)(2*num_values - 1)));
        fmt->sz_accum_rel = 1 + (fmt->pr_accum_rel > 0) + fmt->pr_accum_rel;
      }
    else
      { /* The cumulative counts are half-integers divided by {tot_wt_count}: */
        fmt->pr_accum_rel = (tot_wt_count <= 1 ? 1 : (int32_t)digits(2*tot_wt_count - 1));
        fmt->sz_accum_rel = 1 + (int32_t)(fmt->pr_accum_rel > 0) + fmt->pr_accum_rel;
      }
  }

options_t *get_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem"); 

    if (argparser_keyword_present(pp, "-channel"))
      { if (argparser_keyword_present_next(pp, "R"))
          { o->channel = 0; }
        else if (argparser_keyword_present_next(pp, "red"))
          { o->channel = 0; }
        else if (argparser_keyword_present_next(pp, "G"))
          { o->channel = 1; }
        else if (argparser_keyword_present_next(pp, "grn"))
          { o->channel = 1; }
        else if (argparser_keyword_present_next(pp, "green"))
          { o->channel = 1; }
        else if (argparser_keyword_present_next(pp, "B"))
          { o->channel = 2; }
        else if (argparser_keyword_present_next(pp, "blu"))
          { o->channel = 2; }
        else if (argparser_keyword_present_next(pp, "blue"))
          { o->channel = 2; }
        else 
          { o->channel = (int32_t)argparser_get_next_int(pp, 0, 3); }
      }
    else
      { o->channel = -1; }
      
    if (argparser_keyword_present(pp, "-mask"))
      { o->wt_filename = argparser_get_next_non_keyword(pp); }
    else
      { o->wt_filename = NULL; }
    
    /* Allocate and parse ignored values: */
    o->ignore = bool_vec_new(PNM_FILE_MAX_MAXVAL);
    for (uint32_t v = 0;  v < o->ignore.ne; v++) { o->ignore.e[v] = FALSE; }
    while (argparser_keyword_present(pp, "-ignore"))
      { int32_t v = (int32_t)argparser_get_next_int(pp, 0, PNM_FILE_MAX_MAXVAL);
        if (o->ignore.e[v])
          { fprintf(stderr, "warning: sample value %d already ignored\n", v); }
        else
          { o->ignore.e[v] = TRUE; }
      }
      
    o->verbose = argparser_keyword_present(pp, "-verbose");

    argparser_skip_parsed(pp);
    
    if (argparser_next(pp) != NULL)
      { o->in_filename = argparser_get_next(pp); }
    else
      { o->in_filename = "-"; }
    
    argparser_finish(pp);

    return o;
  }
