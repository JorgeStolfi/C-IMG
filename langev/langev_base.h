/* Basic data types for {langev.c} */
/* Last edited on 2008-01-09 15:14:44 by stolfi */ 

#ifndef langev_base_H
#define langev_base_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <values.h>

#include <bool.h>
#include <vec.h>

#include <langev_lang.h>
#include <langev_gene.h>
#include <langev_name.h>
#include <langev_geog.h>

#define INF INFINITY   

/* STATE OF THE WORLD */

#define langev_base_frames_INFO \
  "  The state of the world at any given instant is described by" \
  " a /frame/.  The main components of a frame are two rectangular" \
  " arrays {L} (/lang/) and {G} (/gene/), with one element per site.  Each element of {L}" \
  " is a 32-bit binary word that models the language" \
  " spoken by the people living at the corresponding site.  Each element" \
  " of {G} is also a 32-bit binary word that models the genetic" \
  " makeup of those people.  A special element value {Ø} is used to mark" \
  " /vacant/ (uninhabited) sites in both arrays.\n" \
  "\n" \
  "  Besides the {L} and {G} arrays, a frame contains other" \
  " time-varying data, which will be detailed further on."
  
#define langev_base_family_INFO \
  "  The simulator assumes that, in each generation, each site" \
  " is linguistically and genetically homogeneous.  Therefore, we" \
  " may as well assume that each site of the world that is not vacant is" \
  " occupied by a single adult /family/, consisting of a couple of twins," \
  " brother and sister, who have the same genetic makeup and who speak" \
  " the same language."

typedef struct sibs_t /* Attributes of a person or twin couple. */
  { lang_t lang;  /* Language, or {NULL_LANG}. */
    gene_t gene;  /* Genome, or {NULL_GENE}. */
    name_t name;  /* Identity (birth certificate number), or {NULL_NAME}. */
    name_t ndad;  /* Identity of father, or {NULL_NAME}. */
    name_t nmom;  /* identity of mother, or {NULL_NAME}. */
  } sibs_t;
  
#define NULL_SIBS \
  (sibs_t){ \
    /*lang:*/ NULL_LANG, \
    /*gene:*/ NULL_GENE, \
    /*name:*/ NULL_NAME, \
    /*ndad:*/ NULL_NAME, \
    /*nmom:*/ NULL_NAME \
  }

typedef struct site_t
  { sibs_t oc;     /* Characteristics of the occupying twins, of {NULL_SIBS} if vacant. */
    double remt;   /* Remoteness of site, used in {langev_move.h}. */
  } site_t; 
  /* An entry of a frame that corresponds to a site. */

typedef struct frame_t /* A simulation frame (world sate). */ 
  { int cols;       /* Number of columns (pixels per row). */
    int rows;       /* Number of rows. */
    site_t *site;   /* Linearized sites: {rows*cols}. */
    double *xtv;    /* X migratory drift per language group (size: N_GR_LANG). */
    double *ytv;    /* Y migratory drift per language group (size: N_GR_LANG) */
    double *prs;    /* Prestige increment per language group (size: N_GR_LANG). */
  } frame_t;

frame_t *new_frame(int nx, int ny);
  /* Allocates a new frame structure, for a world with {nx} columns 
    and {ny} rows of sites. Also allocates the language-related 
    attribute tables.  The sites and tables are left uninitialized. */

site_t *get_site_address(frame_t *frm, site_id_t p);
  /* Returns the address of the site {p} in the frame {frm}.
    The index {x} is implicitly reduced modulo {frm->cols}. */

bool_t do_settle_twins
  ( frame_t *frm,
    geography_t *geo,   /* World's geography. */
    sibs_t *ss,
    site_id_t pp, 
    bool_t verbose
  );
  /* Forcibly places a couple of twin children with attributes {ss}
    at the site {pp} of frame {frm}. If the site is on water, the
    new occupants die, the site remains vacant, and the procedure
    returns FALSE. Otherwise the site gets occupied by the given
    couple, and the procedure returns TRUE. The previous occupants of
    the site in frame {frm}, if any, are murdered in cold blood. */

void fill_sites(frame_t *frm, sibs_t sibs, double remt);
  /* Stores the given [sibs} attributes into all sites of {frm}. */

double get_language_prestige(frame_t *frm, lang_t v);
  /* Returns the numeric prestige of language {v} in state {frm}. */  

void get_language_drift_vector(frame_t *frm, lang_t v, double *xt, double *yt);
  /* Stores in {*xt,*yt} the migration trend vector associated to 
    language {v} in state {frm}. */  

void output_frame
  ( char *prefix, 
    int iep, 
    frame_t *frm, 
    geography_t *geo
  );
  /* Writes a colorized version of {frm} to files called "{prefix}-{NNNNNN}-L.ppm"
    and "{prefix}-{NNNNNN}-G.ppm" where {NNNNNN} is the six-digit decimal
    representation of {iep}. */

/* PRINTOUT */

void print_site(FILE *wr, int nx, int ny, char *pre, sibs_t *sp, site_id_t *pp, double *altp, char *suf);
void print_sibs(FILE *wr, char *pre, sibs_t *sp, char *suf);
void print_coord(FILE *wr, char *pre, int z, char *suf);
void print_site_coords(FILE *wr, int nx, int ny, char *pre, site_id_t p, char *suf);
void print_site_id(FILE *wr, char *pre, site_id_t p, char *suf);
void print_coords(FILE *wr, char *pre, int x, char *sep, int y, char *suf);
void print_altitude(FILE *wr, char *pre, double alt, char *suf);
  /* These procedures print the given info to {wr}. 
    Each pointer field is suppressed if it is NULL.
    
    The {z} parameter is supposed to be either a X or a Y coordinate.
    
    The procedures {print_site} and {print_site_coords} will print the
    coordinates of {p}, rather than its id, assuming grid size {nx} by
    {ny}, using "(" and ")" as coordinate delimiters and " " as the
    separator. For the procedure {print_coords}, the delimiters and
    separator must be part of {pre}, {suf}, and {sep}. */ 

#endif
