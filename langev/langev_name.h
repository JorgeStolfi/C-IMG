/* Tools for simulated genomes */
/* Last edited on 2008-01-06 19:27:47 by stolfi */ 

#ifndef langev_name_H
#define langev_name_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <values.h>

#include <bool.h>
#include <vec.h>

#define N_BITS_NAME 30
  /* Number of traits in each genome. */

typedef uint32_t name_t; 
  /* An identifier for a twin children couple. Values between 0 and
    {2^N_BITS_NAME-1} represent valid names. The two high-order bits are
    used occasionally as marks or to indicate null values. */

#define NULL_NAME (1u << N_BITS_NAME)
  /* An invalid genome code that indicates a vacant site (including water). */
    
#define EDEN_NAME (0u)
  /* The genome of Eden. */

name_t next_name(void);
  /* Returns an original name, which is not {NULL_NAME}. Names may
    repeat after {2^N_BITS_NAME} calls. */

/* PRINTOUT */

void print_name(FILE *wr, char *pre, name_t n, char *suf);
  /* Prints the name {n} to {wr}, in base 26 (or the string 
    "--" if {g == NULL_NAME}). */

#endif
