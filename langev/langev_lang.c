/* See langev_lang.h. */
/* Last edited on 2024-12-21 14:00:26 by stolfi */

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
#include <frgb.h>
#include <frgb_ops.h>

#include <langev_lang.h>
#include <langev_util.h>

double abs_lang_difference(lang_t va, lang_t vb)
  { /* Count differing traits: */
    lang_t mask = 1u;
    double dif = 0;
    int i;
    for (i = 0; i < N_TRAITS_LANG; i++)
      { /* Compare one trait: */
        if ((va & mask) != (vb & mask)) { dif += i+1; }
        /* Go to next trait: */
        mask <<= 1;
      }
    return dif;
  }

double rel_lang_difference(lang_t va, lang_t vb)
 { 
   return abs_lang_difference(va,vb)/MAX_DIFF_LANG;
 }

lang_t mix_languages(lang_t v[], double w[], int n)
  { /* Compute the total weight {wtot}: */
    int k;
    double wtot = 0;
    for (k = 0; k < n; k++) { wtot += w[k]; }
    
    /* Pick each trait at random: */
    lang_t vx = 0u;
    lang_t mask = 1u;
    int i;
    for (i = 0; i < N_TRAITS_LANG; i++)
      { /* Select the donor of trait {i}: */
        int kd = pick_elem(w, n, 0, wtot);
        /* Copy trait {i} from {v[kd]}: */
        vx |= (mask & v[kd]);
        /* Go to next trait: */
        mask <<= 1;
      }
    return vx;
  }

lang_t mutate_language(lang_t v, double prob)
  { /* Start with language {v}: */
    lang_t vm = v;
    /* Flip traits at random: */
    lang_t mask = 1u;
    int i;
    for (i = 0; i < N_TRAITS_LANG; i++)
      { /* Flip trait {i} with probability {prob/(i+1)}: */
        if(drandom() < prob/(i+1)) { vm ^= mask; }
        /* Go to next trait: */
        mask <<= 1;
      }
    return vm;
  }

/* LANGUAGE GROUPS */

#define MINOR_BITS (N_TRAITS_LANG - N_TRAITS_GR_LANG)
#define GROUP_BITS (N_TRAITS_GR_LANG)

#define MINOR_DISP 0
#define GROUP_DISP (MINOR_DISP + MINOR_BITS)

#define MINOR_MASK (((1u << MINOR_BITS) - 1u) << MINOR_DISP)
#define GROUP_MASK (((1u << GROUP_BITS) - 1u) << GROUP_DISP)

int language_group(lang_t v)
  { 
    /* Return the {N_TRAITS_GR_LANG} higher-order traits, as an integer: */
    return (int)((v & GROUP_MASK) >> GROUP_DISP); 
  }

lang_t representative_lang(int k)
  {
    /* Return a language with the given higher-order traits: */
    return (((lang_t)k) << GROUP_DISP) & GROUP_MASK; 
  }

/* COLORING LANGUAGES */

#define HUE_BITS GROUP_BITS
  /* Number of bits in language  that are mapped to hue. */

#define SAT_BITS 2
  /* Number of bits in language that are mapped to saturation. */

#define LUM_BITS (N_TRAITS_LANG - HUE_BITS - SAT_BITS)
  /* Number of bits in language that are mapped to brightness. */

#define MIN_HUE 0.20
#define MAX_HUE 0.80
  /* Min and max brightness. */

#define MIN_LUM 0.40
#define MAX_LUM 0.60
  /* Min and max brightness. */

void map_language_to_color(lang_t v, float rgb[])
  { 
    /* We will scan the traits from highest to lowest: */
    lang_t mask = 1u << (N_TRAITS_LANG - 1);
    /* Select a hue {hue} based on the hightest {HUE_BITS} traits: */
    double p = 0.5; /* Trait positon weight. */
    int j;
    double hue = 0;
    for (j = 0; j < HUE_BITS; j++)
      { if ((v & mask) != 0) { hue += p; }
        /* Reduce {p} by more than half, to reflect trait importance: */
        p = 0.40*p;
        /* Go to next trait: */
        mask >>= 1;
      }
    /* Remap the hue to avoid blues: */
    hue = MIN_HUE + (MAX_HUE - MIN_HUE)*hue;
    /* fprintf(stderr, "    hue = %6.4f\n", hue); */
    
    /* Select a relative saturation {sat} based on the next {SAT_BITS} bits: */
    p = 0.25; /* Trait position weight. */
    double sat = 1.0;
    for (j = 0; j < SAT_BITS; j++)
      { if ((v & mask) != 0) { sat -= p; }
        /* Reduce {p} by more than half, to reflect trait importance: */
        p = 0.40*p;
        /* Go to next trait: */
        mask >>= 1;
      }
        
    /* Apply a hash function on the remaining bits: */
    lang_t vy = v;
    for (j = 0; j < LUM_BITS; j++)
      { if ((vy & 1u) != 0) { vy ^= 19504615u; }
        vy >>= 1;
      }
    mask = (1u << LUM_BITS) - 1; /* Selects the last {} bits of lang. */
    vy = vy & mask;
    /* Select a relative brightness {lum} based on the hashed bits: */
    double lum = MIN_LUM + (MAX_LUM - MIN_LUM)*((double)vy)/((double)mask);
    
    /* Convert {hue,sat,lum} to RGB: */
    frgb_t f = (frgb_t){{ hue, sat, lum}};
    assert((f.c[0] >= 0) && (f.c[0] <= 1.0));
    assert((f.c[1] >= 0) && (f.c[1] <= 1.0));
    assert((f.c[2] >= 0) && (f.c[2] <= 1.0));
    frgb_from_HTY(&f);
    assert((f.c[0] >= 0) && (f.c[0] <= 1.0));
    assert((f.c[1] >= 0) && (f.c[1] <= 1.0));
    assert((f.c[2] >= 0) && (f.c[2] <= 1.0));
    rgb[0] = f.c[0]; 
    rgb[1] = f.c[1]; 
    rgb[2] = f.c[2]; 
  }

#define MAX_WIDTH_LANG 9
  /* Max printed width of a {lang_t} value. */

void print_lang(FILE *wr, char *pre, lang_t v, char *suf)
  {
    if ((pre != NULL) & ((*pre) != 0)) { fputs(pre, wr); }
    if (v == NULL_LANG)
      { int i;
        for (i = 0; i < MAX_WIDTH_LANG; i++) { fputc('-', wr); }
      }
    else
      { /* Print group and details in hex alpha: */
        print_hex_alpha(wr, (v & GROUP_MASK) >> GROUP_DISP, GROUP_BITS); 
        fputc('.', wr);
        print_hex_alpha(wr, (v & MINOR_MASK) >> MINOR_DISP, MINOR_BITS); 
      }
    if ((suf != NULL) & ((*suf) != 0)) { fputs(suf, wr); }
  }
