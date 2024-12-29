#define PROG_NAME "langev"
#define PROG_DESC "simulate language evolution over generations"
#define PROG_VERS "1.0"

#define langev_C_COPYRIGHT "Copyright © 2006 by the State University of Campinas (UNICAMP)"

/* Last edited on 2024-12-21 14:00:42 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

#include <jsfile.h>
#include <jsrandom.h>
#include <argparser.h>
#include <affirm.h>
#include <bool.h>
#include <vec.h>
#include <r2.h>

#include <langev_info.h>
#include <langev_lang.h>
#include <langev_gene.h>
#include <langev_name.h>
#include <langev_base.h>
#include <langev_geog.h>
#include <langev_util.h>
#include <langev_move.h>

typedef struct sim_params_t
  { /* Parameters affecting events of individuals: */
    distr_params_t *brides;  /* Prob distr for bride selection. */
    distr_params_t *tutors;  /* Prob distr of language tutors. */
    distr_params_t *travel;  /* Prob distr for migration before settlement. */
    distr_params_t *settle;  /* Prob distr for settlement site selection. */
    double lang_flip_prob;   /* Prob of random lang trait mutation. */
    double gene_link_prob;   /* Prob of linkage between adjacent gene traits. */
    double gene_flip_prob;   /* Prob of random gene trait mutation. */
    double avg_children;     /* Avg number of children for medium-prestige individuals. */
    int max_settle_trials;   /* Maximum number of settlement attempts per family. */
    /* Parameters affecting the evolution of language attributes: */
    double prs_change;       /* Max prestige change per generation. */
    double mtv_change;       /* Max migration trend change per generation. */
    double mtv_max_len;      /* Max length of migration trend vector. */
  } sim_params_t; 
  /* A collection of parameters that controls the simulation. */

typedef struct options_t /* The command line arguments. */
  { char *relief;    /* Name of input relief file, or NULL. */
    int level;       /* Pixval of sea level. */
    int generations; /* Number of generations to simulate. */
    int step;        /* Output one every this many epochs. */
    int nx;          /* Declared width, or -1 if none.  */
    int ny;          /* Declared height, or -1 if none.  */
    int startx;      /* Column of initial family.  */
    int starty;      /* Row of initial family.  */
    char *prefix;    /* Output file name prefix. */
  } options_t;

typedef struct counts_t /* Counts of things. */
  { int n_men;              /* Men considered for marriage. */
    int n_marriages;        /* Marriages performed. */
    int n_settle_try;       /* Attempts to settle twin couples. */
    int n_settle_suc;       /* Twin couples settled with no deaths. */
    int n_settle_dis;       /* Twin couples displaced after settling. */
    int n_settle_die;       /* Twin couples died while trying to settle. */
    /* The fields {n_cd_XXX} count sites used for {pick_elem}. */
    /* The fields {n_vs_XXX} count sites visited while gethering those {n_cd_XXX} sites. */
    int n_cd_brides, n_vs_brides; /* Potential brides. */                       
    int n_cd_tutors, n_vs_tutors; /* Language tutors. */                        
    int n_cd_settle, n_vs_settle; /* settlement sites. */
    int n_cd_travel, n_vs_travel; /* Migration destinations. */                 
} counts_t;

/* INTERNAL PROTOTYPES */

int main(int argc, char* argv[]);

options_t *parse_options(int argc, char **argv);

void misc_tests(void);

void clear_counts(counts_t *ct);
  /* Clears all fields of {*ct}. */

void print_counts(FILE *wr, counts_t *ct);
  /* Prints the counts {ct} and derived statistics to {stderr}. */

sim_params_t *choose_sim_params(options_t *o);
 /* Returns a record of simulation parameters {sim}, possibly 
   depending on the command line options {o}. */

void update_population_attributes
  ( frame_t *old,
    frame_t *new,
    sim_params_t *spm   /* Simulation parameters. */
  );
  /* Computes the language-related attributes of the new state {new},
    as small random changes from those of the old state {old},
    according to the parameters in {spm}. */
     
void init_language_attributes
  ( sim_params_t *spm,  /* Simulation parameters. */
    int na, 
    double prs[],
    double xtv[], 
    double ytv[]
  );
  /* Initializes {prs[0..na-1]}, {xtv[0..na-1]}, and {ytv[0..na-1]} to random values,
    according to the parameters in {spm}. */

/* 
  SIMULATION PROCEDURES

  All of the {simulate_XXX} procedures below increment the counters in
  the {ct} argument (which must be initialized by the client). */

void simulate_one_generation
  ( frame_t *old,        /* Old frame. */
    int old_iep,         /* Epoch index of old frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    sim_params_t *spm,   /* Simulation parameters. */
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    counts_t *ct         /* Counts of things. */
  );
  /* Given the state of the world {old} at the old epoch {old_iep},
    stores into {new} the state at epoch {old_iep+1}.  */

void simulate_marriages_of_man
  ( frame_t *old,        /* Old frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    sim_params_t *spm,   /* Simulation parameters. */
    site_id_t pf,        /* Site id of father. */
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    counts_t *ct         /* Counts of things. */ 
  );
  /* Marries the man living at site {pf} of {old} with all nearby
    women that the prestige {prs[vf]} of its langauge {vf} will allow.
    
    Then simulates the life of the the pair of twins born of each marriage, 
    from language learning, migration, and settlement in the {new} frame.
    
    The vectors {Q} and {W} are used as work areas; they should be
    allocated by the client, with arbitrary size, and will be expanded
    as needed by the procedure. */

void simulate_life_of_twin_couple
  ( frame_t *old,        /* Old frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    sim_params_t *spm,   /* Simulation parameters. */
    site_id_t pf,        /* Site id of father. */
    site_id_t pm,        /* Site id of mother. */
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    counts_t *ct         /* Counts of things. */
  );
  /* Simulates the life, up to adulthood, of a couple of twin children
    born from the marriage of the man living at site {pf} with
    the woman living at {pm}, both in the {old} frame.
    
    The simulation includes birth, learning the language from their
    parents the local community, pre-marriage migration, and
    settlement in the new frame.
    
    The vectors {Q} and {W} are used as work areas; they should be
    allocated by the client, with arbitrary size, and will be expanded
    as needed by the procedure. */
    
lang_t simulate_language_learning
  ( frame_t *old,
    site_list_t *L,
    double prob_flip
  );
  /* Simulates the language learning process for a couple of twin
    children.
    
    Each trait of language {v} is independently learned from the
    language spoken by some neighbor. At the end, each trait is
    randomyly flipped with probability {prob_flip}.
    
    The neighbors are {qt[0..nw-1]}, and the probability of neighbor
    {qt[k]} being selected as the donor of a trait is proportional
    to{wt[k]}. */

bool_t simulate_migration_and_settlement
  ( frame_t *old,
    frame_t *new,
    geography_t *geo,            /* World's geography. */
    distr_params_t *dpm_settle,  /* Site distribution parameters for settlement. */
    distr_params_t *dpm_travel,  /* Site distribution parameters for migration. */
    sibs_t *sb,                  /* (READONLY) Attributes of twin couple. */
    site_id_t pb,                /* Site id of birth site. */
    int max_trials,              /* Max trials while trying to settle one couple. */
    site_list_t *L,              /* Work area. */
    pqueue_t *Q,                 /* Work area. */
    counts_t *ct                 /* Counts of things. */ 
  );
  /* Simulates the possible migration and attempted settlement of a
    new pair of twin children with attributes {*sb} that were born at
    site {pb}.
    
    The procedure first tries to settle the children at their birth
    site {pb} or thereabouts, by calling {simulate_settlement_attempt}.
    If that attempt fails, or succeeds by displacing a previously
    settled couple, the procedure calls {simulate_migration} to select
    a site elsewhere for the still-unsettled couple, and then tries to
    settle them at that site or thereabouts, again by calling
    {simulate_settlement_attempt}.
    
    After {max_trials} attempts, any children that are still unsettled
    just die.
    
    The initial site must be on dry land, and the children will be
    settled only on dry land. The parameters {dpm_settle} define the
    probability distribution for site selection in
    {simulate_settlement_attempt}. The parameters {dpm_travel} define the
    probability distribution of steps for {simulate_migration}.
    
    The procedure returns TRUE iff the new couple was settled without
    any other couple dying. */

void simulate_migration
  ( frame_t *old,
    frame_t *new,
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Site distribution parameters. */
    double xt,           /* X coord of migration trend vector. */  
    double yt,           /* Y coord of migration trend vector. */
    sibs_t *ss,          /* Attributes of migrating couple (readonly). */
    site_id_t *pp,       /* Site id of initial site (IN) and final site (OUT).  */ 
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    counts_t *ct         /* Counts of things. */ 
  );
  /* Simulates the primary migration of a new pair of twin children
    with attributes {*ss} who start at site {*pp}. The end point
    of the migration is returned in {*pp}.  The children are not
    settled anywhere.
    
    The initial site must be on dry land, and the final site will be
    on dry land. The probability of each final site is affected by the
    distribution parameters {dpm} and the migration trend vector
    {(xt,yt)}. */

void simulate_settlement_attempt
  ( frame_t *old,
    frame_t *new,
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Site distribution parameters for settlement. */
    sibs_t *sp,          /* Attributes of twin couple (IN) and displaced couple (OUT). */
    site_id_t *pp,       /* Site id of initial site (IN) and final site (OUT).  */ 
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    counts_t *ct         /* Counts of things. */
  );
  /* Tries to place a couple of twin children with characteristics
    {*sp} in the frame {new}, at site {*pp} or thereabouts. 
    
    The procedure first selects the candidate settlement site
    {pr} with a distribution centered at {*pp} and defined
    by the attributes {*sp} and the parameters {dpm}. Let {*sr} be the
    pointer to that site's present occupants.
    
    If the selected site is vacant (i.e. {*sr} is a null {sibs_t}), 
    the procedure swaps {*sp} with {*sr}, settling the couple there. 
    
    If the selected site is occupied, the procedure decides whether
    the the couple {*sp} should replace the present occupants {*sr}.
    The probability of that happening depends on the prestiges of the
    languages spoken by {*sp} and by {*sr} in the {new} frame. If the
    {*sp} couple wins, the procedure swaps {*sp} with {*sr};
    otherwise, it leaves {*sp} and {*sr} undisturbed.
    
    Either way, the coordinates of the selected site {pr} are copied
    into {*pp}. So, if {*sp} is not null after the call, it is a
    homeless couple, sitting desolate by the road before the gate of
    site {*pp} */

bool_t simulate_random_settlement
  ( frame_t *old,
    frame_t *new,
    geography_t *geo,            /* World's geography. */
    distr_params_t *dpm_settle,  /* Site distribution parameters for settlement. */
    distr_params_t *dpm_travel,  /* Site distribution parameters for settlement. */
    sibs_t *sb,                  /* (READONLY) Attributes of couples to settle. */
    int max_trials,              /* Max settlement attempts per new couple. */
    site_list_t *L,              /* Work area. */
    pqueue_t *Q,                 /* Work area. */
    counts_t *ct                 /* Counts of things. */
  );
  /* Tries to place a couple of twin children with individual
    attributes {*sb} at a randomly chosen site of frame {new}.
    
    Specifically, picks a dry land random site {pr} in {new},
    pretends that the couple {*sb} was born there, and calls
    {simulate_migration_and_settlement} at that site; with the
    candidate site probability distribuiton defined by {dpm_settle},
    the migration distribution defined by {dpm_travel}, and
    {max_trials} maximum settlement trials.
    
    Returns TRUE iff the new couple was settled without any other
    couple dying. */

/* ROUTINES */
 
static bool_t DEBUG = FALSE;

int main(int argc, char* argv[])
  { 
    /* Parse the command line arguments: */
    options_t *o = parse_options(argc, argv);
    
    /* misc_tests(); */
    
    /* Obtain the world's geography, and get its size {(nx,ny)}: */
    geography_t *geo = new_geography(o->relief, o->level, o->nx, o->ny);
    int nx = geo->nx, ny = geo->ny;
    int nsites = nx*ny;
    
    /* Create the simulation parameter tables: */
    sim_params_t *spm = choose_sim_params(o);
    
    /* Create two consecutive frames {old} and {new}: */
    frame_t *old = new_frame(nx, ny);
    frame_t *new = new_frame(nx, ny);

    /* Initialize the old frame to all empty: */
    fill_sites(old, NULL_SIBS, +INF);
    
    /* Initialize the language-related attributes of {old}: */
    init_language_attributes(spm, N_GR_LANG, old->prs, old->xtv, old->ytv);
    
    /* Settle Adam and his sister Eve in the {old} frame: */
    { sibs_t adev = 
        (sibs_t){
          /*lang:*/ EDEN_LANG,
          /*gene:*/ EDEN_GENE,
          /*name:*/ EDEN_NAME,
          /*ndad:*/ NULL_NAME,
          /*nmom:*/ NULL_NAME 
        };
      fprintf(stderr, "settling the first twins at [%4d %4d]\n", o->startx, o->starty);
      site_id_t pe = site_id(o->startx,o->starty, nx,ny); /* Eden site id. */
      bool_t res = do_settle_twins(old, geo, &adev, pe, FALSE);
      demand(res, "failure to settle the initial couple"); 
    }
    
    /* Work areas for {simulate_one_generation}: */
    site_list_t *L = site_list_new(nsites);
    pqueue_t *Q = pqueue_new();
    pqueue_realloc(Q, /* nmax: */ nsites, /* zlim:*/ nsites);
    
    /* Allocate the counts vector: */
    counts_t ct;
    
    /* Simulate the requested generations: */
    int old_iep = 0; /* Index of epoch of {old} */
    while (TRUE)
      {
        /* Output the old frame: */
        if (old_iep % o->step == 0)
          { output_frame(o->prefix, old_iep, old, geo); }
        
        /* Do we need one more generation? */
        if (old_iep >= o->generations) { break; }
        
        /* Apply the generation procedure: */
        clear_counts(&ct);
        simulate_one_generation(old, old_iep, new, geo, spm, L, Q, &ct);
        
        /* Print statistics: */
        print_counts(stderr, &ct);
        
        /* Recycle the {old} frame, redefine the {new} frame as {old}: */
        { frame_t *tmp = old; old = new; new = tmp; }
        
        old_iep++;
      }
    return(0);
  }

void init_language_attributes
  ( sim_params_t *spm,  /* Simulation parameters. */
    int na, 
    double prs[],
    double xtv[], 
    double ytv[]
  )
  { int k; 
    for (k = 0; k < na; k++)
      { /* Prestige is a random real value in {{-1 _ +1]}: */
        prs[k] = 2*drandom() - 1;
        /* Migratory drift is a random vector in the unit disk of radius {MAX_DRIFT}: */
        r2_t d;
        r2_throw_ball(&d);
        xtv[k] = d.c[0]*spm->mtv_max_len;
        ytv[k] = d.c[1]*spm->mtv_max_len;
      }
  }

void clear_counts(counts_t *ct)
  {
    ct->n_men = 0;              /* Men considered for marriage. */
    ct->n_marriages = 0;        /* Marriages performed. */
    ct->n_settle_try = 0;       /* Attempts to settle twin couples. */
    ct->n_settle_suc = 0;       /* Twin couples settled with no deaths. */
    ct->n_settle_dis = 0;       /* Twin couples displaced after settling. */
    ct->n_settle_die = 0;       /* Twin couples died while trying to settle. */
            
    ct->n_cd_brides = ct->n_vs_brides = 0; /* Potential brides. */
    ct->n_cd_tutors = ct->n_vs_tutors = 0; /* Language tutors. */
    ct->n_cd_settle = ct->n_vs_settle = 0; /* settlement sites. */
    ct->n_cd_travel = ct->n_vs_travel = 0; /* Migration destinations. */
  }

void print_counts(FILE *wr, counts_t *ct)
  {
    fprintf(wr, "total counts:\n");
    fprintf(wr, "\n");
    fprintf(wr, "%10d men considered for marriage.\n", ct->n_men);
    fprintf(wr, "%10d marriages performed.\n", ct->n_marriages); 
    fprintf(wr, "%10d attempts to settle twin couples.\n", ct->n_settle_try); 
    fprintf(wr, "%10d twin couples settled with no deaths.\n", ct->n_settle_suc); 
    fprintf(wr, "%10d twin couples displaced after settling.\n", ct->n_settle_dis);
    fprintf(wr, "%10d twin couples died while trying to settle.\n", ct->n_settle_die);
    fprintf(wr, "\n");
    
    fprintf(wr, "averages per marriage and birth:\n");
    fprintf(wr, "\n");

    double avg_settle_try = ((double)ct->n_settle_try)/((double)ct->n_marriages);
    double avg_settle_suc = ((double)ct->n_settle_suc)/((double)ct->n_marriages);
    double avg_settle_dis = ((double)ct->n_settle_dis)/((double)ct->n_marriages);
    double avg_settle_die = ((double)ct->n_settle_die)/((double)ct->n_marriages);

    fprintf(wr, "%8.2f attempts to settle twin couples.\n", avg_settle_try); 
    fprintf(wr, "%8.2f twin couples settled with no deaths.\n", avg_settle_suc); 
    fprintf(wr, "%8.2f twin couples displaced after settling.\n", avg_settle_dis); 
    fprintf(wr, "%8.2f twin couples died while trying to settle.\n", avg_settle_die); 
    fprintf(wr, "\n");

    fprintf(wr, "\n");
    fprintf(wr, "total candidate sites (and sites visited to get them):\n");
    fprintf(wr, "\n");
    fprintf(wr, "%10d (%10d) potential brides.\n",        ct->n_cd_brides, ct->n_vs_brides); 
    fprintf(wr, "%10d (%10d) language tutors.\n",         ct->n_cd_tutors, ct->n_vs_tutors); 
    fprintf(wr, "%10d (%10d) settlement sites.\n",        ct->n_cd_settle, ct->n_vs_settle); 
    fprintf(wr, "%10d (%10d) migration destinations.\n",  ct->n_cd_travel, ct->n_vs_travel); 
    fprintf(wr, "\n");
    
    fprintf(wr, "average candidate sites per marriage and birth:\n");
    fprintf(wr, "\n");
    
    double avg_cd_brides = ((double)ct->n_cd_brides)/((double)ct->n_marriages);
    double avg_cd_tutors = ((double)ct->n_cd_tutors)/((double)ct->n_marriages);
    double avg_cd_settle = ((double)ct->n_cd_settle)/((double)ct->n_marriages);
    double avg_cd_travel = ((double)ct->n_cd_travel)/((double)ct->n_marriages);
    
    fprintf(wr, "%8.2f potential brides.\n",       avg_cd_brides); 
    fprintf(wr, "%8.2f language tutors.\n",        avg_cd_tutors); 
    fprintf(wr, "%8.2f settlement sites.\n",       avg_cd_settle); 
    fprintf(wr, "%8.2f migration destinations.\n", avg_cd_travel); 
  }

void simulate_one_generation
  ( frame_t *old,
    int old_iep,
    frame_t *new,
    geography_t *geo,    /* World's geography. */
    sim_params_t *spm,   /* Simulation parameters. */
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    counts_t *ct         /* Counts of things. */ 
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    int nsites = nx*ny; 
    
    /* Clear out sibs data of new frame: */
    fill_sites(new, NULL_SIBS, +INF);

    /* Slowly modify population attributes: */
    update_population_attributes(old, new, spm);

    /* Try to marry everyone and process their children: */
    int nmen = 0;
    site_id_t p;
    for (p = 0; p < nsites; p++)
      { site_t *cxy = get_site_address(old, p);
        if (cxy->oc.lang != NULL_LANG) 
          { /* DEBUG = (nmen == old_iep); */
            simulate_marriages_of_man(old, new, geo, spm, p, L, Q, ct);
            nmen++;
          }
      }

    /* Count occupied sites in the new frame: */
    { int n = 0;
      for (p = 0; p < nsites; p++)
        { site_t *cxy = get_site_address(new, p);
          if (cxy->oc.lang != NULL_LANG) { n++; }
        }
      fprintf(stderr, "alive = %d\n", n);
    }
  }

void update_population_attributes
  ( frame_t *old,
    frame_t *new,
    sim_params_t *spm   /* Simulation parameters. */
  )
  { bool_t debug = DEBUG;
    int na = N_GR_LANG;
    int k; 
    for (k = 0; k < na; k++)
      { /* Twiddle the prestige terms: */
        double pk = old->prs[k];
        pk += spm->prs_change * (2*drandom() - 1);
        /* Clip to {[-1 _ +1]}, by folding: */
        if (pk > +1.00) { pk = +2.00 - pk; }
        if (pk < -1.00) { pk = -2.00 - pk; }
        new->prs[k] = pk;

        /* Twiddle the drift vector: */
        r2_t u; 
        r2_throw_ball(&u);
        double dx = old->xtv[k] + spm->mtv_change * u.c[0];
        double dy = old->ytv[k] + spm->mtv_change * u.c[1];
        double m = hypot(dx, dy);
        if (m > spm->mtv_max_len) 
          { double s = spm->mtv_max_len/m; dx *= s; dy *= s; }
        new->xtv[k] = dx;
        new->ytv[k] = dy; 
        
        if (debug)
          { fprintf(stderr, "group %03d", k);
            fprintf(stderr, "  prs = %+7.4f", new->prs[k]);
            fprintf(stderr, "  dvt = %+6.3f %+6.3f", new->xtv[k], new->ytv[k]);
            fprintf(stderr, "\n");
          }
      }
  }
  
void simulate_marriages_of_man
  ( frame_t *old,       /* Old frame. */
    frame_t *new,       /* New frame. */
    geography_t *geo,   /* World's geography. */
    sim_params_t *spm,  /* Simulation parameters. */
    site_id_t pf,       /* Site id of father. */
    site_list_t *L,     /* Work area. */
    pqueue_t *Q,        /* Work area. */
    counts_t *ct        /* Counts of things. */
  )
  { 
    bool_t debug = DEBUG;
    /* bool_t debug = (ct->n_men == 0); */
              
    int nx = geo->nx;
    int ny = geo->ny;
    
    demand(get_altitude(geo, pf) > 0.0, "father is a big fish?");
    
    ct->n_men++;

    site_t *cf = get_site_address(old, pf); /* Home site of father. */
    sibs_t *sf = &(cf->oc);
    assert(sf->lang != NULL_LANG);

    /* Compute the father language's prestige: */
    double prsf = get_language_prestige(old, sf->lang);

    double ach = spm->avg_children;     /* Avg num of children of max-prestige individuals. */
    double nch = ach * exp(M_LN2*prsf); /* Number of children to generate. */

    if (debug) 
      { print_site(stderr, nx, ny, "marrying man ", sf, &pf, NULL, "");
        fprintf(stderr, " prs = %+6.2f nch = %5.2f\n", prsf, nch);
      }

    while ((nch >= 1.0) || (nch > drandom()))
      { 
        /* Select the mother among adult girls near {pf} that
          speak languages similar to the man's.  
          cannot use {select_site_by_walking}: if population density
          is low, a bounded random walk may fail to find a bride. */
        site_id_t pm = pf;
        int nQ, mQ;
        select_site_from_mob_matrix(old, new, geo, spm->brides, sf, &pm, L, &nQ, &mQ);
        ct->n_cd_brides += nQ; ct->n_vs_brides += mQ;
        
        /* Grab the mother's site {cm} and language {vm}: */
        site_t *cm = get_site_address(old, pm);   /* Mother's home site. */
        sibs_t *sm = &(cm->oc); 
        if (debug) { print_site(stderr, nx, ny, "  selected mother ", sm, &pm, NULL, "\n"); }
        
        /* Simulate marriage and birth: */
        simulate_life_of_twin_couple(old, new, geo, spm, pf, pm, L, Q, ct);
        ct->n_marriages++;
        
        if (debug) { fprintf(stderr, "\n"); }
        
        /* One more child sired: */
        nch = nch - 1;
      }
  }
  
void simulate_life_of_twin_couple
  ( frame_t *old,       /* Old frame. */
    frame_t *new,       /* New frame. */
    geography_t *geo,   /* World's geography. */
    sim_params_t *spm,  /* Simulation parameters. */
    site_id_t pf,       /* Site id of father. */
    site_id_t pm,       /* Site id of mother. */
    site_list_t *L,     /* Work area. */
    pqueue_t *Q,        /* Work area. */
    counts_t *ct        /* Counts of things. */
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    
    bool_t debug = DEBUG;
    
    site_t *cf = get_site_address(old, pf);
    site_t *cm = get_site_address(old, pm);

    /* Choose birth site {pb} -- a parent's home, either dad's or mom's: */
    bool_t born_at_dads = (drandom() < 0.5);
    site_id_t pb = (born_at_dads ? pf : pm);

    /* Compute the children's genome {gb}: */
    gene_t gf = cf->oc.gene; /* Father's genome. */
    gene_t gm = cm->oc.gene; /* Mother's genome. */
    gene_t gb = compute_childrens_genome(gf, gm, spm->gene_link_prob, spm->gene_flip_prob);
    
    /* Choose birth language {vb} (either dad's or mom's): */
    lang_t vb = (drandom() < 0.5 ? cf->oc.lang : cm->oc.lang);

    /* Choose the children's name {mc}: */
    name_t mb = next_name();
    /* Assemble the record {sb} of the children's attributes at birth: */
    sibs_t sb = 
      (sibs_t){
        /*lang:*/ vb,
        /*gene:*/ gb,
        /*name:*/ mb,
        /*ndad:*/ cf->oc.name,
        /*nmom:*/ cm->oc.name 
      }; 
    if (debug) { print_site(stderr, nx, ny, "  children born   ", &sb, &pb, NULL, "\n"); }

    /* Gather adults living near {b} who speak languages similar to {vb}: */
    int mQ;
    collect_sites_from_mob_matrix(old, new, geo, spm->tutors, &sb, pb, L, &mQ); 
    ct->n_cd_tutors += L->n; ct->n_vs_tutors += mQ;

    /* Mix the birth language {vb} with the languages of that community: */
    lang_t vc = simulate_language_learning(old, L, spm->lang_flip_prob);

    /* Add some random mutation: */
    vc = mutate_language(vc, spm->lang_flip_prob);
    
    /* Assemble the record {sc} with the children's grown-up attributes: */
    sibs_t sc = 
      (sibs_t){
        /*lang:*/ vc,
        /*gene:*/ gb,
        /*name:*/ mb,
        /*ndad:*/ cf->oc.name,
        /*nmom:*/ cm->oc.name 
      }; 
    if (debug) { print_lang(stderr, "  children learn language  ", sc.lang, "\n"); }
    
    /* Try to settle children: */
    bool_t ok =
      simulate_migration_and_settlement
        ( old, new, geo, spm->settle, spm->travel, &sc, pb, spm->max_settle_trials, L, Q, ct );
    if (debug) { fprintf(stderr, "  population %s\n", (ok ? "increased" : "remained the same")); }
  }

bool_t simulate_migration_and_settlement
  ( frame_t *old,
    frame_t *new,
    geography_t *geo,            /* World's geography. */
    distr_params_t *dpm_settle,  /* Site distribution parameters for settlement. */
    distr_params_t *dpm_travel,  /* Site distribution parameters for migration. */
    sibs_t *sb,                  /* (READONLY) Attributes of twin couple. */
    site_id_t pb,                /* Site id of birth site. */
    int max_trials,              /* Max trials while settling one couple. */
    site_list_t *L,              /* Work area. */
    pqueue_t *Q,                 /* Work area. */
    counts_t *ct                 /* Counts of things. */ 
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    
    bool_t debug = (pb == site_id(310,80, nx,ny));
    /* bool_t debug = DEBUG; */ 
    
    /* The current couple is {sc}: */
    sibs_t sc = (*sb); 
    demand(sc.lang != NULL_LANG, "cannot settle nobody");
    
    /* The candidate sites for settlement are initially those around {pb}: */
    site_id_t pc = pb;
    
    /* Perform {max_trials} attempts at settlement and migration: */
    int trials = 0; /* Number of settlement attempts: */
    while (trials < max_trials)
      { /* Make one attempt to settle the {sc} couple near {pc}: */
        if (debug) { print_site(stderr, nx, ny, "  trying to settle", &sc, &pc, NULL, "\n"); }
        simulate_settlement_attempt(old, new, geo, dpm_settle, &sc, &pc, L, Q, ct);
        trials++;
        ct->n_settle_try ++;
        
        if (sc.lang == NULL_LANG)
          { /* Settlement succeeded completely: */ 
            if (debug) { print_site_coords(stderr, nx, ny, "  settlement succeded at ", pc, "\n"); }
            ct->n_settle_suc ++;
            return TRUE;
          }
        else
          { /* Settlement failed, or displaced a previous couple: */ 
            if (debug) { print_site(stderr, nx, ny, "  leaves unsettled", &sc, &pc, NULL, "\n"); } 
          }
        
        /* Compute the migration trend vector {(xt,yt)} for the couple {sc}: */
        double xt, yt;
        get_language_drift_vector(new, sc.lang, &xt, &yt);

        /* Try to migrate from {pc}, more or less in the direction {(xt,yt)}: */
        simulate_migration(old, new, geo, dpm_travel, xt, yt, &sc, &pc, L, Q, ct);
        if (debug) { print_site(stderr, nx, ny, "  migrate children", &sc, &pc, NULL, "\n"); }
      }
    
    /* Too many attempts, must give up: */
    ct->n_settle_die ++;
    if (debug) { print_site(stderr, nx, ny, "  discard children", &sc, &pc, NULL, "\n"); }
    return FALSE;
  }

void simulate_migration
  ( frame_t *old,
    frame_t *new,
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Site distribution parameters. */
    double xt,           /* X coord of migration trend vector. */  
    double yt,           /* Y coord of migration trend vector. */
    sibs_t *ss,          /* Attributes of migrating couple (readonly). */
    site_id_t *pp,       /* Site id of initial site (IN) and final site (OUT).  */ 
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    counts_t *ct         /* Counts of things. */ 
  )
  { 
    int nQ, mQ;
    select_site_by_walking(old, new, geo, dpm, xt, yt, ss, pp, L, &nQ, &mQ);
    ct->n_cd_travel += nQ; ct->n_vs_travel += mQ;
  }

void simulate_settlement_attempt
  ( frame_t *old,
    frame_t *new,
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Site distribution parameters for settlement. */
    sibs_t *sp,          /* Attributes of twin couple (IN) and displaced couple (OUT). */
    site_id_t *pp,       /* Site id of initial site (IN) and final site (OUT).  */ 
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    counts_t *ct         /* Counts of things. */
  )
  { 
    int nx = geo->nx;
    int ny = geo->ny;
    
    bool_t debug = ((*pp) == site_id(310,80, nx,ny));
    /* bool_t debug = DEBUG; */
    
    if (debug) { print_site(stderr, nx, ny, "    try to settle ", sp, pp, NULL, "\n"); }

    /* Select the settlement site, put that into {*pp}: */
    int nQ, mQ;
    select_site_from_mob_matrix(old, new, geo, dpm, sp, pp, L, &nQ, &mQ);
    ct->n_cd_settle += nQ; ct->n_vs_settle += mQ;

    /* Get pointer {*sr} to its present occupants: */
    site_t *cr = get_site_address(new, *pp); 
    sibs_t *sr = &(cr->oc);
    if (debug) { print_site(stderr, nx, ny, "    selected site ", sp, pp, NULL, "\n"); }

    /* Decide which couple keeps the site: */
    bool_t swap; 
    if (sr->lang == NULL_LANG)
      { /* Site is vacant: */
        swap = TRUE;
      }
    else
      { /* Site is occupied -- compare prestiges: */
        double prsp = get_language_prestige(new, sp->lang);
        double prsr = get_language_prestige(new, sr->lang);
        double dp = (prsp - prsr)/2; /* Rel diff in prestiges. */
        assert(fabs(dp) <= 1.0);
        double prwin = (1 + dp*dp*dp)/2;
        /* Compute outcome of dispute: */
        swap = drandom() < prwin;
        ct->n_settle_dis ++;
      }
    if (debug) { fprintf(stderr, "    new couple %s\n", (swap ? "settles" : "bounces")); }
      
    /* If the newly arrived couple wins, swap with current occupants: */
    if (swap) 
      { /* Swap the children {sc} and the present occupants: */
        sibs_t st = *sp; *sp = *sr; *sr = st;
      }
  }

lang_t simulate_language_learning
  ( frame_t *frm,
    site_list_t *L,
    double prob_flip
  )
  { 
    int n = L->n; /* Count of tutors. */
    /* Gather tutor languages {vt[0..nw-1]} from {frm}: */
    lang_t vt[n];
    int k;
    for (k = 0; k < n; k++)
      { site_id_t pk = L->itm[k];
        site_t *sk = get_site_address(frm, pk);
        vt[k] = sk->oc.lang;
      }
      
    /* Compute the mean language of the community: */
    lang_t vc = mix_languages(vt, L->val, n);
    
    /* Mutate the language: */
    vc = mutate_language(vc, prob_flip);
    
    return vc;
  }
  
bool_t simulate_random_settlement
  ( frame_t *old,
    frame_t *new,
    geography_t *geo,            /* World's geography. */
    distr_params_t *dpm_settle,  /* Site distribution parameters for settlement. */
    distr_params_t *dpm_travel,  /* Site distribution parameters for settlement. */
    sibs_t *sb,                  /* (READONLY) Attributes of couples to settle. */
    int max_trials,              /* Max settlement attempts per new couple. */
    site_list_t *L,              /* Work area. */
    pqueue_t *Q,                 /* Work area. */
    counts_t *ct                 /* Counts of things. */
  )
  { 
    int nx = geo->nx;
    int ny = geo->ny;

    site_id_t pb;
    do {
      int xb = int32_abrandom(0, nx);
      int yb = int32_abrandom(0, ny);
      pb = site_id(xb,yb, nx,ny);
    } while (get_altitude(geo, pb) <= 0.0);
    
    return
      simulate_migration_and_settlement
        ( old, new, geo, dpm_settle, dpm_travel, sb, pb, max_trials, L, Q, ct );
  }

#define MAX_LEVEL (uint16_image_MAX_SAMPLE - 1)
  /* The maximum sea level, so that one can have some dry land. */

#define MAX_GENERATIONS 100000000
  /* The maximum number of simulated generations (safety). */

options_t *parse_options(int argc, char **argv)
  { 
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t*)notnull(malloc(sizeof(options_t)), "out of mem");
    
    /* Parse keyword parameters: */
    
    argparser_get_keyword(pp, "-generations");
    o->generations = argparser_get_next_int(pp, 0, MAX_GENERATIONS);
    
    if (argparser_keyword_present(pp, "-step"))
      { o->step = argparser_get_next_int(pp, 1, MAX_GENERATIONS); }
    
    o->nx = -1;
    o->ny = -1;
    o->relief = NULL;
    o->level = 127;
    if (argparser_keyword_present(pp, "-size"))
      { o->nx = argparser_get_next_int(pp, 1, MAX_SIZE_GEOG);
        o->ny = argparser_get_next_int(pp, 1, MAX_SIZE_GEOG);
      }
    
    if (argparser_keyword_present(pp, "-relief"))
      { o->relief = argparser_get_next(pp);
        o->level = argparser_get_next_int(pp, 0, MAX_LEVEL);
      }
    
    if ((o->relief == NULL) && (o->nx == -1))
      { argparser_error(pp, "must specify at least one of \"-size\" or \"-relief\""); }
    
    argparser_get_keyword(pp, "-start");
    o->startx = argparser_get_next_int(pp, 0, MAX_SIZE_GEOG-1);
    o->starty = argparser_get_next_int(pp, 0, MAX_SIZE_GEOG-1);
    
    argparser_get_keyword(pp, "-prefix");
    o->prefix = argparser_get_next(pp);

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Check for spurious arguments: */
    argparser_finish(pp);

    return o;
  }

sim_params_t *choose_sim_params(options_t *o)
  {
    /* Allocate the command line argument record: */
    sim_params_t *spm = (sim_params_t*)notnull(malloc(sizeof(sim_params_t)), "out of mem");
    
    /* Parameters affecting the evolution of language attributes: */
    spm->prs_change =  0.10;  /* Max prestige change per generation. */
    spm->mtv_change =  0.20;  /* Max migration trend change per generation. */
    spm->mtv_max_len = 3.00;  /* Max length of migration trend vector. */
    
    /* Parameters affecting events of individuals: */
    spm->lang_flip_prob = 0.02; /* Prob of random lang trait mutation. */
    spm->gene_link_prob = 0.75; /* Prob of linkage between adjacent gene traits. */
    spm->gene_flip_prob = 0.02; /* Prob of random gene trait mutation. */
    spm->avg_children = 1.25;   /* Avg number of children for medium-prestige individuals. */
    spm->max_settle_trials = 3; /* Maximum number of settlement attemps. */
    
    /* Distribution for bride selection: */
    spm->brides = new_distr_params
      ( 
        /*shw_rel_diff:*/    0.5,    /* Bride search over shallow water is OK. */
        /*dpw_rel_diff:*/    +INF,   /* Bride search over deep water is impossible. */
        /*mean_trip_remt:*/  2.0,    /* Search brides over about this many cells. */
        /*start_site_wt:*/   0.001,  /* Avoid marriages with man's own sister. */
        /*use_new_state:*/   FALSE,  /* Brides are chosen from adult population. */
        /*mean_lang_diff:*/  0.10,   /* Tolerable relative lang diff of bride. */
        /*vacant_factor:*/   0.00    /* Cannot marry non-existing woman. */
      );
      
    /* Distribution of childhood language tutors: */
    spm->tutors = new_distr_params
      ( 
        /*shw_rel_diff:*/    3.0,    /* Hard to chat across shallow water. */
        /*dpw_rel_diff:*/    +INF,   /* Impossible to chat across deep water. */
        /*mean_trip_remt:*/  1.0,    /* Tutors this many cells away are half as influent. */
        /*start_site_wt:*/   1.000,  /* Parents are important tutors. */
        /*use_new_state:*/   FALSE,  /* Tutors come from the adult population. */
        /*mean_lang_diff:*/  0.30,   /* Tolerable relative lang diff of tutor. */
        /*vacant_factor:*/   0.00    /* Cannot learn much from non-existing tutors. */
      );
      
    /* Distribution of candidate sites for settlement: */
    spm->settle = new_distr_params
      ( 
        /*shw_rel_diff:*/    2.0,   /* Searching for sites across shallow water is so-so. */
        /*dpw_rel_diff:*/    +INF,  /* Searching for sites across deep water is impossible. */
        /*mean_trip_remt:*/  1.0,   /* Search for sites by about this many cells. */
        /*start_site_wt:*/   1.000, /* Birthplace is a good site to settle in. */
        /*use_new_state:*/   TRUE,  /* Settlement site selection depends on new state. */
        /*mean_lang_diff:*/  INF,   /* Language of settlement dest. site is irrelevant. */
        /*vacant_factor:*/   1.00   /* But vacant sites are quite attractive. */
      );
      
     /* Distribution for target site of migration after failed settlement attempt: */
    spm->travel = new_distr_params
      ( 
        /*shw_rel_diff:*/    0.25,  /* Migration across across shallow water is easy. */
        /*dpw_rel_diff:*/    20.0,  /* Migration across across deep water is nearly impossible. */
        /*mean_trip_remt:*/  1.0,   /* Randomly drift by this many cells. */
        /*start_site_wt:*/   0.010, /* No use retrying where we were bounced/displaced. */
        /*use_new_state:*/   TRUE,  /* Migration dest. selection depends on new state. */
        /*mean_lang_diff:*/  INF,   /* Language of migration dest. site is irrelevent. */
        /*vacant_factor:*/   1.00   /* Vacant sites are quite attractive. */
      );
      
   return spm;
  }

void misc_tests(void)
  { int ix;
    for (ix = 0; ix <= 100; ix++)
      { double x = 0.1*ix;
        fprintf(stderr, "x = %4.1f", x);
        fprintf(stderr, "  exp = %24.22f", exp(-x));
        fprintf(stderr, "  erf = %24.22f", erf(x));
        fprintf(stderr, "  truncated_exp = %24.22f", truncated_exp(-x));
        fprintf(stderr, "  truncated_sigmoid = %24.22f", truncated_sigmoid(x));
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
  }

