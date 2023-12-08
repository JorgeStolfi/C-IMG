/* See langev_name.h. */
/* Last edited on 2008-01-07 22:37:03 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

#include <jsrandom.h>
#include <affirm.h>
#include <bool.h>
#include <vec.h>

#include <langev_name.h>

static name_t last_name = EDEN_NAME;

name_t next_name(void)
  {
    last_name = (last_name + 1) % (1 << N_BITS_NAME);
    return last_name;
  }

#define MAX_WIDTH_NAME 7
  /* Max length of a {name_t} in base 26. */

void print_name(FILE *wr, char *pre, name_t n, char *suf)
  {
    if ((pre != NULL) & ((*pre) != 0)) { fputs(pre, wr); }
    if (n == NULL_NAME)
      { int i;
        for (i = 0; i < MAX_WIDTH_NAME; i++) { fputc('-', wr); }
      }
    else
      { /* Print the name in base 26, low-end-first: */
        /* Print the first letter in uppercase: */
        fputc('A' + (n % 26), wr); n /= 26;
        /* Print the remaining letters in lowercase, supressing trailing 'a's: */
        int k;
        for (k = 1; k < MAX_WIDTH_NAME; k++) 
          { if (n == 0)
              { fputc(' ', wr); }
            else
              { fputc('a' + (n % 26), wr); }
            n /= 26;
          }
      }
    if ((suf != NULL) & ((*suf) != 0)) { fputs(suf, wr); }
  }
