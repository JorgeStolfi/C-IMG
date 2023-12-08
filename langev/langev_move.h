/* Tools for selecting sites with various distributions */
/* Last edited on 2008-01-11 02:29:06 by stolfi */ 

#ifndef langev_move_H
#define langev_move_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <values.h>

#include <bool.h>
#include <vec.h>
#include <pqueue.h>
#include <double_spmat.h>

#include <langev_lang.h>
#include <langev_gene.h>
#include <langev_name.h>
#include <langev_base.h>
#include <langev_geog.h>
#include <langev_list.h>

/* TRAVEL DIFFICULTY */
  
#define langev_move_travel_INFO \
  "  Several sub-processes of the generation procedure (such as" \
  " marriage, language learning, and settlement) involve travel between" \
  " sites.  Travel always follows an 8-connected path, that is, a sequence" \
  " of chess-king moves.\n" \
  "\n" \
  "  The likelyhood of any travel is strongly affected by the" \
  " /difficulty/ of the trip; which depends not only on its length" \
  " but also on the altitudes and slopes encountered along the path.\n" \
  "\n" \
  "  The simulator aims to model the world before reliable oceanic" \
  " navigation, so trips are largely confined to dry land or shallow" \
  " water.  Trips over deep water are assumed to be exceedingly unlikely, and" \
  " are generally ignored by the simulator.\n" \
  "\n" \
  "  Travel over shallow water may be considered easier or more" \
  " difficult than travel on dry land, depending on the nature" \
  " of the trip, as detailed in the following sections.  "
  
double step_difficulty
  ( geography_t *geo,  /* World's geography. */
    site_id_t pa,      /* Site id of origin site. */
    int xs,            /* X increment of step. */
    int ys,            /* Y increment of step. */
    double shw_fac,    /* Relative difficulty of travel on shallow water. */
    double dpw_fac     /* Relative difficulty of travel to deep water. */
  );
  /* The difficulty of the step from site {pa} to the site reached from
    {pa} by the step {(xs,ys)}.
    
    The step {(xs,ys)} must be a king move, and the origin and
    destination sites must exist; otherwise the result is {+INF}. The
    origin coordinates {pa} must have been reduced to the range
    {{0..nx-1} × {0..ny-1}}.
    
    The parameters {shw_fac} and {dpw_fac} are factors that
    multiplying the difficulty of steps in shallow water and deep
    water, respectively. See {langev_move_stdiff_INFO} for details. */
    
#define langev_move_stdiff_INFO \
  "  The difficulty of an elementary step depends basically" \
  " on its length and on the altitudes of the two sites.  A step" \
  " of unit length on dry flat terrain, just above sea level, has" \
  " difficulty 1.  The difficulty increases with increasing altitude" \
  " and with increasing slope (upward or downward).\n" \
  "\n" \
  "  If the initial and final site of the step are both at zero" \
  " altitude (shallow water), the difficulty of the step gets" \
  " multiplied by a /shallow water factor/.  If both" \
  " sites have negative altitude (deep water), the difficulty" \
  " gets multiplied by a /deep water factor/.  Depending on" \
  " the type and purpose of the trip, travel on water" \
  " may be easier (factor {< 1}) or harder (factor {> 1}) than" \
  " on dry land.  Usually, however, travel in deep water is more" \
  " difficult than in shallow water.\n" \
  "\n" \
  "  If *only one* of the sites is on deep water, of any depth, " \
  " its altitude is assumed to be zero (shallow water) for" \
  " the purposes of computing its difficulty.  That is, a" \
  " step that moves from dry land or shallow water into deep" \
  " water (or vice-versa) is no harder than one that moves" \
  " into shallow water (or vice-versa).  The deep water factor" \
  " applies only if *both* sites are in deep water."

/* RANDOM SITE CHOICES 
  
  The simulator often needs to randomly choose a site for some
  purpose (such as marriage, migration, settlement, etc.).
  
  The probability of choosing a site {pk} is usually a function of
  several variables:
  
    * the altitude of {pk}, in particular whether it is 
    on dry land or in water;
    
    * the people that are living at {pk} in some epoch, in particular
    whether it is vacant or occupied;
   
    * the attributes {ss} of the /selector/, who is the person or twin
    couple who is supposedly choosing the site;
    
    * the remoteness of {pk} from some /initial site/ {pi}, e.g. the 
    current abode of the selector;
     
    * the length and direction of a /migration trend vector/, that
    depends on the selector's language and varies over time.
    
  The probability distribution also depends on the purpose of the trip
  (marriage, language acquisition, learning, etc).  Thus, for each
  purpose there is a different set of /distribution parameters/ which
  determine the precise effect of those variables on the
  probabilities.

  PROBABILITY WEIGHT FUNCTIONS
  
  The probability of choosing a site is computed by collecting a
  relatively small set of /candidate sites/, assigning a non-negative
  /weight/ to each site, and normalizing those weights to unit sum.
  
  The weight {weight(pk)} of a site {pk} is computed as the product of
  two quantities, the /remoteness weight factor/ {remt_wtf(pk)} and the
  /site quality weight factor/ {site_wtf(pk)}.  Both factors are always 
  real numbers in the range {[0 _ 1]}. 
  
  The remoteness factor {remt_wtf(pk)} depends exclusively on the
  difficulty of reaching {pk} from the initial site {pi}. The function
  {remt_wtf} must be monotonic non-increasing, and fall very quickly
  to 0 as {rk} increases. The candidate sites {cands(pi)} are
  precisely those for which {remt_wtf} is non-zero.
  
  The site quality factor function {site_wtf(pk)} may depend on
  arbitrary data; including, for example, the attributes of people who
  are or were settled at {pk}. This factor need not decrease with
  remoteness, and often varies subtantially from site to site. In
  particular, sites on water generally have {site_wtf(pk) == 0}.
  However, site_wtf(pk) must be positive for at least one site in
  {cands(pi)}
  
  Now suppose that the candidate sites are listed in order of
  increasing remoteness, and {pz} is the last site enumerated. Then we
  know that the weight of any candidate {pw} yet to be examined will
  satisfy {weight(pw) <= remt_wtf(pz)}. Thus, we can stop the
  enumeration once the sum of {weight(pk)} for the candidates already
  listed in  ??? Improve and fix.
  
*/

typedef struct distr_params_t
  { /* Parameters that affect the computuation of the difficulty of trips: */
    double shw_rel_diff;   /* Remoteness factor for trips over shallow water. */
    double dpw_rel_diff;   /* Remoteness factor for trips over deep water. */
    /* Parameters that affect the remoteness weight factor function {remt_wtf}: */
    double mean_trip_remt; /* Destination remoteness that reduces the weight by 1/2. */
    /* Parameters that affect the site quality weight factor function {site_wtf}: */
    double start_site_wt;  /* Weight when destination site is equal to the starting site. */
    bool_t use_new_state;  /* FALSE looks at old state, TRUE at new state. */
    double mean_lang_diff; /* Rel. language difference that reduces the weight by 1/2. */
    double vacant_factor;  /* Weight multiplier for vacant sites. */
    /* Precomputed mobility matrices: */
    mob_matrix_t *D;       /* Elementary steps and their difficulties. */
    mob_matrix_t *M;       /* Reachable sites and their remotenesses. */
  } distr_params_t;
  /* The paramter {mean_trip_remt} is the remoteness {rm} such that
    {remt_wtf(rm)} is {1/2}.  Likewise, the parameter {mean_lang_diff} 
    is the relative language difference {rv} that causes {site_wtf(pk)}
    to be halved, all other factors being equal.
    
    The parameters {mean_lang_diff} and {vacant_factor} 
    may refer to the status of the site in the {old} or {new} state, depending on 
    {use_new_state}. */
    

distr_params_t *new_distr_params
  ( double shw_rel_diff, 
    double dpw_rel_diff, 
    double mean_trip_remt, 
    double start_site_wt,     
    bool_t use_new_state,      
    double mean_lang_diff,
    double vacant_factor
  );
  /* Returns a pointer {dpm} to a newly allocated distributon parameters record,
    with given field values.   The mobility matrices {dpm->D} and {dpm->M} are
    set to NULL. */

/* 
*/

double compute_site_remoteness_weight_factor
  ( double rk,           /* Its remoteness from the base site.  */ 
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
  );
  /* Computes  probability weight {wk} for site {pk}, whose remoteness is {rk}.
    The probability distribution is defined by the parameters {dpm,ss,pi}. */ 

double compute_site_quality_weight_factor
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    site_id_t pk,        /* Id of site in question.  */ 
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t pi         /* Id of starting site.  */ 
  );
  /* Computes the (unnormalized) probability weight {wk} for site {pk}, whose remoteness is {rk}.
    The probability distribution is defined by the parameters {dpm,ss,pi}. */ 

double compute_site_weight
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    site_id_t pk,        /* Id of site in question.  */ 
    double rk,           /* Its remoteness from the base site.  */ 
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t pi         /* Id of starting site.  */ 
  );
  /* Computes the (unnormalized) probability weight {wk} for site {pk}, whose remoteness is {rk}.
    The probability distribution is defined by the parameters {dpm,ss,pi}. */ 

void precompute_step_difficulties(geography_t *geo, distr_params_t *dpm);
  /*  Sets {dpm->D} to the elementary step difficulty matrix for the
    geography {geo}, computed using the motion difficulty parameters
    in {*dpm}. See {compute_step_difficulty_matrix} for details. */

void precompute_reachabilities(geography_t *geo, distr_params_t *dpm);
  /*  Sets {dpm->M} to a reachability matrix for the geography {geo},
    with some cutoff {dmax} derived from the parameters in {dpm}.  The
    cutoff is such that sites with remoteness higher than {dmax} will
    have extremely low weight.  
    
    The elementary steps and their difficulties are taken from the
    elementary step difficulty matrix {dpm->D}. The procedure calls
    {precompute_step_difficulties} if necessary to make {dpm->D != NULL}. */

/* MOBILITY MATRICES 
  
  A {mob_matrix_t} is a /mobility matrix/, a matrix whose rows and
  colums correspond to a certain subset {H} of the sites, and whose
  entries are {double}s (incluing {±INF}) associated to pairs of sites.
  If {M} is a mobility matrix, we denote by 
  {M[A,B]} the value of the element in row {A} and column {B}
  of the martix.

  Mobility matrices are intended for applications where most entries
  are /trivial/ (0 or {+INF}, depending on their nature). Therefore, a
  mobility matrix are stored using sparse-matrix techniques.
  
  A mobility matrix can be viewed also as a real-labeled directed
  graph, whose nodes are a subset {H} of the sites, and whose edges
  are a small subset of {H × H}. In that case, the tail (origin) of
  the edge is the row site, and the head (destination) is the column
  site. Thus {M[A,B]} is the label on the edge that goes from site
  {A} to site {B}.

  Since mobility matrices are usually big and expensive to compute,
  they are generally used for properties that depend only on the
  geography or other rarely-changing data, and not on the state of the
  population. For instance, the mobility matrix may be used to store
  the difficulty of travel between sites. */

typedef struct mob_matrix_t
  { double_spmat_t mat; /* Nontrivial entries, sorted by {row} then {col}. */
    int_vec_t start;    /* Start of entries of each site. */
  } mob_matrix_t;
  /* The {col} and {row} fields in the entries of {mat} are site numbers.
    
    These entries of {mat} are sorted by increasing {row}.  Within each row,
    the order is arbitrary.
    
    The edges coming out of site number {p} are {mat.e[jini..jlim-1]}
    where {jini = start.e[p]} and {jlim = start.e[p+1]}. Thus, the
    vector {start} has {nx*ny+1} entries. In particular, site number
    {p} has no non-trivial outgoing edges if and only if {jini ==
    jlim}. Also, the meaningful entries of {mat} are
    {mat.e[0..nmat-1} where {nmat = start.e[nx*ny]}. Beware that
    {mat->ne} is usually greater than {nmat}. */

mob_matrix_t *mob_matrix_new(int nx, int ny);
  /* Returns a pointer {M} to a newly allocated {mob_matrix_t} record,
    suitable for a world with {nx} columns and {ny} rows.
    
    The integer vector {M->start} is fully allocated, with {nx*ny+1}
    elements. The sparse matrix {M->mat} is allocated assuming 5
    edges per site on average. Thus, when appending entries to
    {M->mat}, be sure to use {double_spmat_expand} on it. The
    contents of {M->start} and {M->mat} are not initialized. */

/* STEP DIFFICULTY MATRICES
  
  In an /elementary/ mobility matrix, the edges are limited to a small
  set of /elementary steps/ between adjacent sites --- either the four
  unit axial moves, of the eight king moves.
  
  A special case of elementary matrix is a /step difficulty/ matrix,
  where each entry {M[A,B]} is the difficulty of the step from {A} to
  {B}. The value is {+INF} if the step is invalid, not elementary, or
  impossible to take; and is 0 if and only if {A == B}. */

mob_matrix_t *compute_step_difficulty_matrix
  ( geography_t *geo,  /* World's geography. */
    double shw_fac,    /* Relative difficulty of travel on shallow water. */
    double dpw_fac     /* Relative difficulty of travel to deep water. */
  );
  /* Computes the step difficulty matrix for the given geography
    and the given water difficulty factors. The entries in each row
    are sorted in order of increasing difficulty. */

/* REACHABILITY MATRICES

  A /reachability matrix with cutoff {dmax}/ is a mobility matrix which
  has an edge from site {A} to {B} if and only if both sites are on
  dry land, and there is a sequence of zero or more elementary steps
  (possibly passing through water sites) whose total difficulty is
  {dmax} or less.  The value of the entry {[A,B]} is the 
  minimum difficulty of any such path.  Edges that are missing 
  are supposed to have value {+INF}.
  
  Such a matrix is useful, for example, when a man living at a site
  {A} has to select a bride.  We assume that he first makes a list of
  all adult women that he can reach from {A} with trips up to a
  certain distance {dmax}. Then he will marry a woman from this set,
  living at {B}, with an arbitrary probability distribution --- which
  may depend on many factors, constant or state-dependent, 
  including perhaps the difficulty of reaching {B} from {A}.
  
  This model may not be very realistic, but has the merit of
  restricting the homes of potential brides to a small and *static*
  subset of all sites. */
  
mob_matrix_t *compute_reachability_matrix
  ( geography_t *geo,  /* World's geography. */
    mob_matrix_t *D,   /* Step difficulty matrix. */
    double dmax       /* Cutoff difficulty. */
  );
  /* Computes the reachability matrix with cutoff {dmax} for the 
    given geography, given the elementary step difficulty matrix {D}.
    
    The procedure is faster if the entries of {D}, within each row,
    are sorted by increasing difficulty. */
  
void gather_reachable_sites
  ( int nx,                     /* X size of site grid. */
    int ny,                     /* Y size of site grid. */
    int pi,                     /* Number of starting site. */
    mob_matrix_t *D,            /* Step difficulty matrix. */
    double dmax,                /* Maximum path difficulty allowed. */
    site_list_t *L,             /* (OUT) list of reachable sites. */
    pqueue_t *Q,                /* (WORK) priority queue. */
    int *nv                     /* (OUT) number of sites visited. */
  );
  /* Finds the sites that are reachable from site number {pi} by a
    sequence of zero or more elementary steps, whose total difficulty
    is at most {dmax}. The valid steps and their difficulties are
    taken from the mobility matrix {D}. Uses Dijkstra's algorithm.
    
    The procedure returns its result in the list {L}, which must be
    allocated by the client. The items of {L} are the reachable sites;
    the value associalted to an item {pk} in {L} is the remoteness of
    {pk} from {pi}.
    
    The client must also preallocate queue {Q}, which will be expanded
    as needed.  
    
    The procedure also returns in {*nv} the count of visited sites, 
    that is, sites that were examined but were beyond the {dmax} 
    threshold. */

/* GATHERING SITES */

typedef double step_diff_func_t(site_id_t pc, int xs, int ys); 
  /* Type of a function that returns the difficulyt of a step
    from site {pc} to a neighboring site {pc+(xs,ys)}.
    
    The result must be strictly positive. If the result is {INF}, that
    step is never used.
    
    For best results, if the step is on level terrain just above sea
    level, the result should be the Euclidean length of the step (1 or
    {sqrt(2)}). */

typedef double remt_weight_func_t(double rc); 
  /* Type of a function that returns a numeric weight for a site of
    remoteness {rc} (that can be reached only by trips with difficulty
    {rc} or more). The result must be non-negative, and must be
    non-increasing with {rc}. For best results, it should tend quickly
    to 0 as {rc} increases. */

typedef double site_weight_func_t(double rc, site_id_t pc); 
  /* Type of a function that returns a numeric weight for the site
    {pc}, that has remoteness {rc}. The result must lie between 0
    and 1. This function will be multplied by a {remt_weight_func_t}
    to give the actual site weight. */

void collect_sites_from_mob_matrix
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t pi,        /* Id of starting site.  */ 
    site_list_t *L,      /* (OUT) list of reachable sites. */
    int *mQ              /* (OUT) number of sites visited. */
  );
  /* Gathers all sites that are reachable from {pi} accoding to the
    reachability matrix {dpm->M}, and computes their probablity
    weights according to a distribution defined by the parameters
    {dpm,ss,pi}, assuming a null migration trend vector {(xt,yt)}.
    
    If the matrix {dpm->M} is NULL, the procedure calls
    {precompute_reachabilities} to build it.
    
    The sites with nonzero weight are returned in {L}, sorted in order
    of increasing remoteness of travel from {pi}. The values
    associated to those items {L} are the corresponding probability
    weights. The weights are non-negative but not normalized.
    
    The procedure also returns in {*mQ} the number of sites that were
    examined while collecting those sites. */

void collect_sites_dijkstra
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    double xt,           /* X coordinate of intended displacement. */
    double yt,           /* Y coordinate of intended displacement. */
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t pi,        /* Id of starting site.  */ 
    site_list_t *L,      /* (OUT) list of reachable sites. */
    pqueue_t *Q,         /* (WORK) priority queue. */
    int *mQ              /* (OUT) number of sites visited. */
  );
  /* Gathers all the sites with relevant weight from a distribution
    defined by the central site {pi} and the parameters {dpm,ss,pi},
    and the trend vector {(xt,yt)}.
    
    The sites found are returned in the list {L}, in order of
    increasing remoteness of travel from {pi}. The associated values
    are the corresponding weights. The weights are non-negative but
    not normalized.
    
    The procedure also returns in {*mQ} the number of sites that were
    visited while collecting those sites. */

void select_site_from_mob_matrix
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t *pp,       /* (IN) Id of reference site; (OUT) Id of selected site. */ 
    site_list_t *L,      /* Work area. */
    int *nQ,             /* (OUT) number of sites considered. */
    int *mQ              /* (OUT) number of sites visited. */
  );
  /* Selects a site randomly among the sites that are reachable from
    {*pp} according to the mobility matrix {dpm->M}. The choice uses a
    probability distribution defined by the parameters {dpm,ss,pp},
    assuming a null trend vector {(xt,yt)}. The selected site is
    returned in {*pp} itself.
    
    If the matrix {dpm->M} is NULL, the procedure calls
    {precompute_reachabilities} to build it.
    
    The procedure also returns in {*nQ} the number of sites in the
    population from which the site was selected, and in {*mQ} the
    number of sites that were visited while collecting those sites. */

void select_site_dijkstra
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    double xt,           /* X coordinate of intended displacement. */
    double yt,           /* Y coordinate of intended displacement. */
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t *pp,       /* (IN) Id of reference site; (OUT) Id of selected site. */ 
    site_list_t *L,      /* Work area. */
    pqueue_t *Q,         /* Work area. */
    int *nQ,             /* (OUT) number of sites considered. */
    int *mQ              /* (OUT) number of sites visited. */
  );
  /* Selects a site randomly from a distribution defined by the 
    parameters {dpm,ss,pp},  and the trend vector {(xt,yt)}.
    The selected site is returned in {*pp} itself.
    
    The procedure also returns in {*nQ} the number of sites in the
    population from which the site was selected; and in {*mQ} the
    number of sites that were visited while collecting those {*nQ}
    sites. */

void select_site_by_walking
  ( frame_t *old,        /* Old frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    double xt,           /* X coordinate of intended displacement. */
    double yt,           /* Y coordinate of intended displacement. */
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t *pp,       /* (IN) Id of reference site; (OUT) Id of selected site. */ 
    site_list_t *L,      /* Work area. */
    int *nQ,             /* (OUT) number of sites considered. */
    int *mQ              /* (OUT) number of sites visited. */
  );
  /* Selects a site by performing a random walk from site {*pp},
    governed by the distribution parameters {dpm,sx} and the
    trend vector {(xt,yt)}.
    
    The procedure also returns in {*nQ} the number of sites in the
    population from which the site was selected; and in {*mQ} the
    number of sites that were visited while collecting those {*nQ}
    sites. */

/* DIJKSTRA'S ALGORITHM */

void generic_dijkstra
  ( geography_t *geo,           /* World's geography. */
    site_id_t pi,               /* Id of starting site. */
    step_diff_func_t step_df,   /* Defines the difficulty of a step. */
    remt_weight_func_t remt_wt, /* Defines the remoteness-related weight of a site. */
    site_weight_func_t site_wt, /* Defines the endsite-related weight factor. */
    site_list_t *L,             /* (OUT) list of reachable sites. */
    pqueue_t *Q,                /* (WORK) priority queue. */
    int *mQ                     /* (OUT) number of sites visited. */
  );
  /* Gathers a set of sites in {frm} that are near {pi} and most 
    relevant in a client-defined sense.
    
    The sites found are returned in the list {L}. The value associated
    to a site is its probability weight. The items are sorted in order
    of increasing remoteness from {pi}. The weights are non-negative
    but not normalized
    
    All the sites returned by the procedure are reachable from
    {pi} by paths consisting entirely of king moves which
    have finite total difficulty. The procedure calls the
    client-provided function {step_df} to compute the difficulty of
    each step; which may be {+INF}, meaning ``do not take this step''.
    In the current implementation, the step {ps} will always be a
    nonzero king move (that is, each increment is {-1}, {0}, or {+1}).
    
    Let {PATH[k]} denote the least-difficult trip from the starting
    site {pi} to a site {pk}. The difficulty of that path is, by
    definition, the remoteness of {pk}.
    
    The weight of a site {pc} is defined as the product of a factor
    {remt_wt(rc)} related to the site's remoteness {rc}, and a factor
    {site_wt(rc,pc)} that may depend on {rc} and on other features of
    the site.
    
    Both functions must return non-negative numbers. The function
    {remt_wt} must be non-increasing with {rc}. The result of {site_wt}
    must lie between 0 and 1.
    
    The procedure stops the enumeration when there are no more
    reachable sites, or when the total weight of the sites still to be
    enumerated (as estimated from {remt_wt}) is negligible compared to
    the total weight of the sites already gathered.
    
    The procedure also returns in {*mQ} the number of sites that were
    visited while collecting the sites. */

/* RANDOM WALK SIMULATOR */

typedef double halt_weight_func_t(double dc, site_id_t pc); 
  /* Type of a function that returns a probability weight for halting
    the walk the site {pc}, that was reached after a walk with
    total difficulty {dc}.
    
    The result must lie in {[0_1]}, but need not be normalized in any
    way. The function may return {0} to mean ``do not halt at this site''. */

typedef double step_weight_funt_t(double es, double ds, site_id_t pc, int xs, int ys); 
  /* Type of a function that returns a probability weight for taking 
    the step from site {pc} to site {pc+(xs,ys)}; given the
    Euclidean length {es} of the step {(xs,ys)}, and its difficulty {ds}. 
    
    The result must lie in {[0_1]}, but need not be normalized in any way. The
    function may return {0} to mean ``do not take this step''; and it
    should do so when the destinantion site does not exist. */

void random_walk
  ( geography_t *geo,           /* World's geography. */
    site_id_t *pp,              /* (IN) Id of starting site; (OUT) Id of final site. */
    step_diff_func_t step_df,   /* Defines the difficulty of a step. */
    halt_weight_func_t halt_wt, /* Defines the probability weight of halting at a site. */
    step_weight_funt_t step_wt, /* Defines the probability weight of taking a step. */
    site_list_t *L,             /* Work area. */
    int *nQ,                    /* (OUT) Number of steps taken along the walk. */
    int *mQ                     /* (OUT) number of sites looked at along the walk. */
  );
  /* Simulates a random walk from site {*pp}, governed by the
    step difficulty function {step_df} and the probability weight
    functions {halt_wt,step_wt},
    
    The path will consist entirely of valid sites, conneted by king
    moves which have finite total difficulty. The procedure calls the
    client-provided function {step_df} to compute the difficulty of
    each step; which may be {+INF}, meaning ``do not take this step''.
    In the current implementation, the step {(xs,ys)} will always be a
    nonzero king move (that is, each increment is {-1}, {0}, or {+1}).
    
    The total difficulty of a walk is, by definition, the sum
    of the difficulties of its steps. 
    
    The procedure calls {halt_wt} and {step_wt} to compute the
    relative probability weight of each alternative (halting or taking
    a step). Both functions must return non-negative numbers;
    {halt_wt} should return 0 if the current site is unacceptable, and
    {step_wt} must return 0 if the proposed step is impossible. These
    weights need not be normalized; the procedure will normalize them
    to obtain the probabilities of the steps.
    
    The procedure stops the enumeration only when the halting alternative
    gets seleted, or when all alternatives have zero weight.  
    To avoid loops, the distribution defined by {dpm} should 
    have {dpm->mean_trip_remt < +INF}.
    
    The procedure also returns in {*nQ} the number of sites in the
    walk (at least 1); and in {*mQ} the number of sites that were
    considered when selecting the next move (which should be 9 times
    {*nQ}, for king moves + halt). */

/* PRINTOUT */

void print_site_list
  ( FILE *wr, 
    frame_t *frm, 
    geography_t *geo, 
    char *hd, 
    site_list_t *L,
    char *ft 
  );
  /* Prints the sites in {Q} and their values (assumed to be
    probability weights) one per line, preceded by the {hd} string and
    terminated by the {ft} string. The procedure also prints the
    current occupant of each site in the frame {frm}. */

void print_pqueue
  ( FILE *wr, 
    frame_t *frm, 
    geography_t *geo, 
    char *hd, 
    pqueue_t *Q,
    char *ft 
  );
  /* Prints the sites in {Q} and their values (assumed to be
    remotenesses) one per line, preceded by the {hd} string and
    terminated by the {ft} string. The procedure also prints the
    current occupant of each site in the frame {frm}. */
    
#endif
