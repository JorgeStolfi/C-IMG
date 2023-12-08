/* See langev_gene.h. */
/* Last edited on 2008-01-07 22:28:52 by stolfi */

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

#include <langev_gene.h>
#include <langev_util.h>

#define FATHER_BITS 8
#define MOTHER_BITS 8
#define MIXING_BITS (N_TRAITS_GENE - FATHER_BITS - MOTHER_BITS)

#define FATHER_DISP (MOTHER_DISP + MOTHER_BITS)
#define MOTHER_DISP (MIXING_DISP + MIXING_BITS)
#define MIXING_DISP 0

#define FATHER_MASK (((1u << FATHER_BITS) - 1u) << FATHER_DISP)
#define MOTHER_MASK (((1u << MOTHER_BITS) - 1u) << MOTHER_DISP)
#define MIXING_MASK (((1u << MIXING_BITS) - 1u) << MIXING_DISP)

gene_t compute_childrens_genome(gene_t gf, gene_t gm, double prlink, double prmut)
  { gene_t gx = mix_genomes(gf, gm, prlink);    /* Mix the mixing parts. */
    gene_t gc = combine_gene_parts(gf, gm, gx); /* Child's genome (pristine). */
    return mutate_genome(gc, prmut);
  }

gene_t father_gene_traits(gene_t g)
  { return (g & FATHER_MASK); }

gene_t mother_gene_traits(gene_t g)
  { return (g & MOTHER_MASK); }

gene_t mixing_gene_traits(gene_t g)
  { return (g & MIXING_MASK); }

gene_t combine_gene_parts(gene_t gf, gene_t gm, gene_t gx)
  {  return 
       (gf & FATHER_MASK) |
       (gm & MOTHER_MASK) |
       (gx & MIXING_MASK);
  }

gene_t mix_genomes(gene_t ga, gene_t gb, double prlink)
  { /* Start with either {ga} or {gb}, equal probs: */
    if (drandom() < 0.5) { gene_t gt = ga; ga = gb; gb = gt; }
    /* Assemble mixed genome {gx}: */
    gene_t gx = 0u;
    gene_t mask = (1u << MIXING_DISP);
    int i;
    for (i = 0; i < MIXING_BITS; i++)
      { /* Copy mixing trait {MIXING_DISP + i} from {ga}: */
        gx |= (mask & ga);
        /* Randomly switch parents: */
        if (drandom() > prlink) { gene_t gt = ga; ga = gb; gb = gt; }
        /* Go to next trait: */
        mask <<= 1;
      }
    return gx;
  }

gene_t mutate_genome(gene_t g, double prmut)
  { /* Turn {prmut} into the probability per trait: */
    prmut /= N_TRAITS_GENE;
    /* Mutate traits: */
    gene_t mask = 1u;
    int i;
    for (i = 0; i < N_TRAITS_GENE; i++)
      { if (drandom() < prmut) { g ^= mask; }
        /* Go to next trait: */
        mask <<= 1;
      }
    return g;
  }

#define MAX_WIDTH_GENE 8
  /* Max printed width of a {gene_t} value. */

void print_gene(FILE *wr, char *pre, gene_t g, char *suf)
  {
    if ((pre != NULL) & ((*pre) != 0)) { fputs(pre, wr); }
    if (g == NULL_GENE)
      { int i;
        for (i = 0; i < MAX_WIDTH_GENE; i++) { fputc('-', wr); }
      }
    else
      { /* Print the three parts in hex alpha: */
        print_hex_alpha(wr, (g & FATHER_MASK) >> FATHER_DISP, FATHER_BITS); 
        fputc('.', wr);
        print_hex_alpha(wr, (g & MOTHER_MASK) >> MOTHER_DISP, MOTHER_BITS); 
        fputc('.', wr);
        print_hex_alpha(wr, (g & MIXING_MASK) >> MIXING_DISP, MIXING_BITS); 
      }
    if ((suf != NULL) & ((*suf) != 0)) { fputs(suf, wr); }
  }

/* void map_genome_to_color(gene_t g, double rgb[]) */
