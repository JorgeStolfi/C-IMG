#define PROG_NAME "ppvip"
#define PROG_DESC "the PPV Image Processor for manipulation of PPV images"
#define PROG_VERS "1.0"

/* Copyright © 2006 by the State University of Campinas (UNICAMP). */
/* See the copyright, authorship, and warranty notice at end of file. */
/* Last edited on 2021-07-08 15:44:10 by jstolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [-seed NUM] \\\n" \
  "    [VAR=]EXPR \\\n" \
  "    ... \n" \
  "  where EXPR is one of:\n" \
  "    NUM\n" \
  "    VAR\n" \
  "    read([FMT:]FNAME)\n" \
  "    write(EXPR,[FMT:]FNAME)\n" \
  "    add(EXPR,EXPR)\n" \
  "    crop(EXPR,AXIS,NUM,NUM)\n" \
  "  VAR is a variable identifier,\n" \
  "  NUM is a numeric constant,\n" \
  "  AXIS is a nonnegative integer,\n" \
  "  FMT is a file format tag (PPM,PPV,...),\n" \
  "  FNAME is a file name, possibly with escaped chars."

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  ???\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in mar/2006 by Jorge Stolfi, IC-UNICAMP."

/* Must define _GNU_SOURCE in order to get {asprintf} */
#define _GNU_SOURCE
#include <math.h>
#include <values.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
  
#include <bool.h>
#include <sign.h>
#include <affirm.h>

#include <ppv_array.h>
#include <ppv_array_read.h>
#include <ppv_array_write.h>

typedef struct options_t 
  { int seed;            /* Seed for random number generators. */
    command_t *cmds;     /* List of commands. */
  } options_t;

int main (int argc, char **argv);
options_t *get_options(int argc, char **argv);
void data_error(int line, char *msg);

int main (int argc, char **argv)
  { 
    options_t *o = get_options(argc, argv);
    
    return 0;
  }

options_t *get_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    
    argparser_finish(pp);

    return o;
  }

void data_error(int line, char *msg)
  {
    fprintf(stderr, "%s:%d: **%s\n", "-", line, msg);
    exit(1);
  }
