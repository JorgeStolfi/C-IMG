/* Tools for simulated genomes */
/* Last edited on 2008-01-06 19:26:49 by stolfi */ 

#ifndef langev_gene_H
#define langev_gene_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <values.h>

#include <bool.h>
#include <vec.h>

#define N_TRAITS_GENE 30
  /* Number of traits in each genome. */

typedef uint32_t gene_t; 
  /* A genome code. Values between 0 and {2^N_TRAITS_GENE-1} represent human
    genomes; each bit is a genetic trait. The two high-order bits
    are used occasionally as marks or to indicate null values. */

#define NULL_GENE (1u << N_TRAITS_GENE)
  /* An invalid genome code that indicates a vacant site (including water). */
    
#define EDEN_GENE (0u)
  /* The genome of Eden. */

/* MAIN SECTIONS OF THE GENOME */

gene_t father_gene_traits(gene_t g);
  /* Returns the part of the genome {g} that is always inherited 
    from the father, with the traits in their original positions. */

gene_t mother_gene_traits(gene_t g);
  /* Returns the part of the genome {g} that is always inherited 
    from the mother, with the traits in their original positions. */

gene_t mixing_gene_traits(gene_t g);
  /* Returns the part of the genome {g} that is randomly inherited 
    from the mother or the father, with the traits in their
    original positions. */

/* SEXUAL RECOMBINATION */

gene_t compute_childrens_genome(gene_t gf, gene_t gm, double prlink, double prmut);
  /* Computes the genome of a new twin couple given the 
    father's genome {gf} and the mother's genome {gm}. 
    The mixing part of the genome is mixed with linkage probability
    {prlink}. */
    
/* LOW-LEVEL GENOME FUNCTIONS */

gene_t combine_gene_parts(gene_t gf, gene_t gm, gene_t gx);
  /* Returns the genome that has the father's part taken from {gf},
     the mother's part taken from {gm}, and the mixing part 
     taken from {gx}. */

gene_t mix_genomes(gene_t ga, gene_t gb, double prlink);
  /* Returns a random mixture of the mixing parts of the genomes {ga}
    and {gb}. The traits are in their normal positions. Non-mixing
    traits are set to zero.
    
    Th mixing traits are copied sequentially; each trait is copied
    from the same parent as the previous trait with probability
    {prlink}, and from the other parent with probability
    {1-prlink}. */

gene_t mutate_genome(gene_t g, double prmut);
  /* Changes each trait of {g} with probability {prmut/N_TRAITS_GENE}. */

void map_genome_to_color(gene_t g, double rgb[]);

/* PRINTOUT */

void print_gene(FILE *wr, char *pre, gene_t g, char *suf);
  /* Prints the genome {g} to {wr}, as a bit vector (or the string 
    "--" if {g == NULL_GENE}). */

#endif
