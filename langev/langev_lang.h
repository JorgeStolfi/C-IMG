/* Tools for simulated languages */
/* Last edited on 2008-01-06 19:21:07 by stolfi */ 

#ifndef langev_lang_H
#define langev_lang_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <values.h>

#include <bool.h>
#include <vec.h>

#define N_TRAITS_LANG 30
  /* Number of traits in each language. */

typedef uint32_t lang_t; 
  /* A language code. Values between 0 and {2^N_TRAITS_LANG - 1} represent
    human languages; each bit is a language trait. The two high-order
    bits are used occasionally as marks or to indicate null values. */

#define NULL_LANG (1u << N_TRAITS_LANG)
  /* An invalid language code that indicates a vacant site (including water). */
    
#define EDEN_LANG (0u)
  /* The code of the language of Eden. */

#define M_H_LANG 3.99499
  /* The sum of {1/k} for {k=1..N_TRAITS_LANG}. */

lang_t mix_languages(lang_t v[], double w[], int n);
  /* Returns a language whose trait number {i} was
    copied from language {v[k]} with probability 
    proportional to {w[k]}. */

lang_t mutate_language(lang_t v, double prob);
  /* Flips trait number {i} of {g} with probability {prob/(i+1)}. Thus the
    expected number of flipped traits will be {prob} times {M_H_LANG}. */

double abs_lang_difference(lang_t va, lang_t vb);
  /*  The difference between the two languages {va,vb}, defined as the
    sum of {i+1} over all {i} such that the two languages differ on
    trait {i} .
    
    If the languages are regarded as bit vectors, this is also the sum
    of {(i+1)*(va[i]-vb[i])^2} , or the sum of {(i+1)*abs(va[i]-vb[i])}, over
    all {i}. */

#define MAX_DIFF_LANG (0.5*(N_TRAITS_LANG*(N_TRAITS_LANG+1)))
  /* The maximum value of {abs_lang_difference(va,vb)} for two
    languages {va,vb}. */

#define AVG_DIFF_LANG (0.5*MAX_DIFF_LANG)
  /* The expected value of {abs_lang_difference(va,vb)} for two random
    languages {va,vb}. */

double rel_lang_difference(lang_t va, lang_t vb);
  /* The relative distance between two languages {va,vb}; defined as
    {abs_lang_difference(va,vb)} divided by {MAX_DIFF_LANG}. The
    maximum result is 1.0, meaning that the languages differ in every
    trait. The expected result for two random languages is 0.5 . */

/* LANGUAGE GROUPS */

#define N_TRAITS_GR_LANG 8
   /* Number of traits that comprise the major language group. */

#define N_GR_LANG (1 << N_TRAITS_GR_LANG)
  /* Number of entries in language-related attribute tables. */

int language_group(lang_t v);
  /* Returns an index in {0..N_GR_LANG-1} depending on the language {v},
    such that languages with the same  are fairly similar to 
    each other according to {abs_lang_difference}. */
    
lang_t representative_lang(int k);
  /* Returns a language {v} that is a canonical representative of
     language group {k}; namely, such that {language_group(v) = k},
     and it is more or less central in that group. */

void map_language_to_color(lang_t v, float rgb[]);
  /* Stores in {rgb[0..2]} a color respresenting the language {v}. 
     
     The most important language traits are mapped to hue differences.
     Secondary traits are mapped to saturation differences. The
     remaining traits are mapped to brightness differences. Many
     languages will map to the same color. */

/* PRINTOUT */

void print_lang(FILE *wr, char *pre, lang_t v, char *suf);
  /* Prints the language {v} to {wr}, as a bit vector (or the string 
    "--" if {v == NULL_LANG}). */

#endif
