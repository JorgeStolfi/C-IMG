/* See {langev_move.h} */
/* Last edited on 2021-07-17 23:28:16 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>

#include <jsrandom.h>
#include <affirm.h>
#include <bool.h>
#include <sign.h>
#include <double_spmat.h>
#include <pqueue.h>
#include <vec.h>

#include <langev_lang.h>
#include <langev_gene.h>
#include <langev_name.h>
#include <langev_base.h>
#include <langev_move.h>
#include <langev_geog.h>
#include <langev_util.h>

double step_difficulty
  ( geography_t *geo,  /* World's geography. */
    site_id_t pa,      /* Site id of origin site. */
    int xs,            /* X increment of step. */
    int ys,            /* Y increment of step. */
    double shw_fac,    /* Relative difficulty of travel on shallow water. */
    double dpw_fac     /* Relative difficulty of travel to deep water. */
  )
  { /* Null steps have zero difficulty: */
    if ((xs == 0) && (ys == 0)) { return 0.0; }
    /* Presently, {step_difficulty} is defined only for king-steps: */
    assert((abs(xs) <= 1) && (abs(ys) <= 1));
    /* Get the world's dimensions: */
    int nx = geo->nx;
    int ny = geo->ny;
    /* Check origin site validity: */
    int xa = site_X_from_id(pa, nx, ny);
    int ya = site_Y_from_id(pa, nx, ny);
    demand((xa >= 0) && (xa < nx), "bad {xa}"); 
    demand((ya >= 0) && (ya < ny), "bad {ya}"); 
    /* Compute the coordinates {(xb,yb)} of the step's destination: */
    int xb = xa + xs;
    int yb = ya + ys;
    /* Reduce the destination to {{0..xn-1} × {0..ny-1}}: */
    while (xb < 0) { xb += nx; }
    while (xb >= nx) { xb -= nx; }
    /* Check existence of the destination site: */
    assert((xb >= 0) && (xb < nx)); /* Guaranteed by reduction. */
    if((yb < 0) || (yb >= ny)) { /* Dest site does not exist: */ return +INF; }
    /* get the id {pb} of the destination site: */
    site_id_t pb = site_id(xb, yb, nx, ny);
    /* Get altitudes at both sites: */
    double ha = get_altitude(geo, pa);
    double hb = get_altitude(geo, pb);
    /* Ignore water depth if one site is wet and the other site is dry or shallow: */
    if ((hb < 0) && (ha >= 0)) { hb = 0.0; }
    if ((ha < 0) && (hb >= 0)) { ha = 0.0; }
    /* Compute the water-related difficulty factors {ud_wat_a,ud_wat_b} for each site: */
    double wat_fac_a = (ha < 0 ? dpw_fac : (ha == 0 ? shw_fac : 1));
    double wat_fac_b = (hb < 0 ? dpw_fac : (hb == 0 ? shw_fac : 1));
    /* Get Euclidean length {es} of displacement: */
    double es = hypot(xs, ys);
    /* Compute the slope-related difficulty {ud_slo} for unit displacement: */
    double dh = hb - ha; /* Altitude increment. */
    double sh = (0.75*fabs(dh) + 0.25*dh)/es; /* Effective slope of terrain. */
    double ud_slo = 4*sh*sh;
    /* Compute the mean altitude-related difficulty {ud_alt} for unit displacement: */
    double ud_alt_a = 4*ha*ha * wat_fac_a;
    double ud_alt_b = 4*hb*hb * wat_fac_b;
    double ud_alt = (ud_alt_a + ud_alt_b)/2;
    /* Compute total difficulty {ud} for unit displacement: */
    double ud = 1.0 + ud_slo + ud_alt;
    /* Multiply by the step's Euclidean length: */
    return es*ud;
  }

mob_matrix_t *mob_matrix_new(int nx, int ny)
  { demand((nx > 0) && (nx <= MAX_SIZE_GEOG), "invalid {nx}");
    demand((ny > 0) && (ny <= MAX_SIZE_GEOG), "invalid {ny}");
    int nsites = nx*ny; /* It should not overflow. */
    mob_matrix_t *M = (mob_matrix_t *)notnull(malloc(sizeof(mob_matrix_t)), "out of mem");
    M->start = int_vec_new(nsites + 1);
    M->mat = double_spmat_new(5*nsites);
    return M;
  }

mob_matrix_t *compute_step_difficulty_matrix
  ( geography_t *geo,  /* World's geography. */
    double shw_fac,    /* Relative difficulty of travel on shallow water. */
    double dpw_fac     /* Relative difficulty of travel to deep water. */
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    int nsites = nx*ny; /* It should not overflow. */
    
    fprintf(stderr, "precomputing a step difficulty matrix\n");
    fprintf(stderr, "  relative travel difficulty in shallow water = %10.6f\n", shw_fac);
    fprintf(stderr, "  relative travel difficulty in deep water =    %10.6f\n", dpw_fac);
    fprintf(stderr, "  %10d sites\n", nsites);
    
    /* Allocate the matrix: */
    mob_matrix_t *D = mob_matrix_new(nx, ny);

    /* Scan the origin sites in order of site number: */
    int pa; /* Current origin site. */
    int mD = 0; /* Number of edges (non-trivial elems) seen so far: */
    D->start.e[0] = 0;  /* Start of edges of site 0: */
    for (pa = 0; pa < nsites; pa++)
      { 
        /* Get the site coordinates: */
        int xa = site_X_from_id(pa, nx,ny);
        int ya = site_Y_from_id(pa, nx,ny);
        
        /* Remember where site {pa} begins in {D->mat}: */
        int jini = mD;
        
        /* Enumerate the elementary steps out of {pa}: */
        int xs, ys; /* Step displacements. */
        for (ys = -1; ys <= +1; ys++)
          { for (xs = -1; xs <= +1; xs++)
              { /* Compute the destination site coords {(xb,yb)}: */
                int xb = xa + xs;
                int yb = ya + ys;
                /* Reduce the destination to {{0..xn-1} × {0..ny-1}}: */
                while (xb < 0) { xb += nx; }
                while (xb >= nx) { xb -= nx; }
                /* Compute the step's difficulty {ds}: */
                double ds = step_difficulty(geo, pa, xs, ys, shw_fac, dpw_fac);
                if (ds < +INF)
                  { /* Difficulty is non-trivial, append to {D->mat}: */
                    double_spmat_expand(&(D->mat), mD);
                    double_spmat_entry_t *eab = &(D->mat.e[mD]);
                    eab->row = pa; 
                    eab->col = site_id(xb, yb, nx, ny);
                    eab->val = ds;
                    mD++;
                  }
              }
          }
        /* Note where the edges of site {pa} end in {D->mat}: */
        int jlim = mD;
        D->start.e[pa+1] = jlim;
        
        /* Sort the entries {D->mat[jini..jlim-1]}: */
        
        auto sign_t cmp_val(double_spmat_entry_t *ea, double_spmat_entry_t *eb); 
          /* Compares entries {*ea} and {*eb} by value. */
        
        double_spmat_sort_entries(D->mat.e, jini, jlim, cmp_val);
        
        sign_t cmp_val(double_spmat_entry_t *ea, double_spmat_entry_t *eb)
          { if (ea->val < eb->val)
              { return -1; }
            else if (ea->val > eb->val)
              { return +1; }
            else
              { return 00; }
          }
      }

    fprintf(stderr, "  %10d non-trivial edges\n", mD);
    fprintf(stderr, "  %10.2f edges per site\n", ((double)mD)/nsites);
    fprintf(stderr, "\n");

    return D;
  }

mob_matrix_t *compute_reachability_matrix
  ( geography_t *geo,  /* World's geography. */
    mob_matrix_t *D,   /* Step difficulty matrix. */
    double dmax        /* Cutoff difficulty. */
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    int nsites = nx*ny; /* It should not overflow. */
    
    fprintf(stderr, "precomputing a reachability matrix\n");
    fprintf(stderr, "  max difficulty of travel = %10.6f\n", dmax);
    fprintf(stderr, "  %10d sites\n", nsites);
    
    /* Allocate the matrix: */
    mob_matrix_t *M = mob_matrix_new(nx, ny);

    /* Create a work queue of pairs {(p,v)}: */
    assert(nsites-1 <= pqueue_ITEM_MAX); 
    pqueue_t *Q = pqueue_new();
    pqueue_realloc(Q, /*nmax:*/ 1000, /*zmax:*/ nsites-1); 
    
    /* Create a list of reachable sites and their remotenesses: */
    site_list_t *L = site_list_new(nsites);
    
    /* Scan the origin sites in order of site number: */
    int pa; /* Current origin site. */
    int mM = 0; /* Number of edges (non-trivial elems) collected: */
    int ndry = 0; /* Number of dry land sites. */
    M->start.e[0] = 0;  /* Start of edges of site 0: */
    for (pa = 0; pa < nsites; pa++)
      { 
        /* Get the site's altitude: */
        double alta = get_altitude(geo, pa);
        
        if (alta > 0.0)
          { /* Origin site is on dry land: */
            ndry++;

            /* Compute the sites reachable from {pa} by Dijkstra's algorithm: */
            int nv;  /* Count of sites visited. */
            gather_reachable_sites(nx, ny, pa, D, dmax, L, Q, &nv);
            int nt = L->n;  /* Count of reachable sites. */
            
            /* Append the reachable sites to {mM}: */
            int k;
            bool_t reaches_self = FALSE; /* Set to TRUE if {pa} reaches {pa}. */
            for (k = 0; k < nt; k++)
              { /* Get the next reachable site {pb}: */
                site_id_t pb = L->itm[k];
                double rb = L->val[k];
            
                /* Get its altitude: */
                double altb = get_altitude(geo, pb);
                
                if ((rb < +INF) && (altb > 0.0))
                  { /* Destination site {pb} is on dry land. */
                    /* Note self-reaching: */
                    if (pb == pa) { reaches_self = TRUE; }
                    /* Append {pb} to {M}: */
                    double_spmat_expand(&(M->mat), mM);
                    double_spmat_entry_t *eab = &(M->mat.e[mM]);
                    eab->row = pa; 
                    eab->col = pb;
                    eab->val = rb;
                    mM++;
                  }
              }

            /* Check for self-reachability: */
            assert(reaches_self); 
          }

        /* Note where the edges of site {pa} end in {M->mat}: */
        int jlim = mM;
        M->start.e[pa+1] = jlim;
      }

    demand(ndry > 0, "world is all water");

    pqueue_free(Q);
    site_list_free(L);
    
    fprintf(stderr, "  %10d non-trivial edges\n", mM);
    fprintf(stderr, "  %10.2f edges per site\n", ((double)mM)/nsites);
    fprintf(stderr, "  %10.2f edges per dry site\n", ((double)mM)/ndry);
    fprintf(stderr, "\n");

    return M;
  }

void gather_reachable_sites
  ( int nx,                     /* X size of site grid. */
    int ny,                     /* Y size of site grid. */
    int pi,                     /* Number of starting site. */
    mob_matrix_t *D,            /* Step difficulty matrix. */
    double dmax,                /* Maximum path difficulty allowed. */
    site_list_t *L,             /* (OUT) list of reachable sites. */
    pqueue_t *Q,                /* (WORK) priority queue. */
    int *nv                     /* (OUT) number of sites visited. */
  )
  { 
    /* bool_t debug = FALSE; */
    
    demand(pi < nx*ny, "bad {pi}");

    /* Flush the queue, and request sorting by increasing {value}: */
    pqueue_reset(Q);         
    pqueue_set_order(Q, +1);

    /* Flush the gathered site list: */
    site_list_reset(L);
    
    /* Initialize the state. */
    (*nv) = 0; /* Number sites visited so far. */

    /* Sites {pt[0..nt-1]} are the /finalized/ sites, namely those for
      which we have already determined the minimum-difficulty path.
      The remoteness of site {pt[i]} (the difficuty of the easiest
      path) is {rt[i]}, for {i} in {0..(*nt)-1}. These sites are sorted in
      increasing order of {rt[i]}.
      
      A site {pf} has been finalized if and only if {kt[pf]} is the index
      of {pf} in the {pt} list; that is, if and only if {kt[pf]} is in {0..(*nt)-1},
      and {pt[kt[pf]] == pf}.  If {pf} has not been finalized, {kt[pf]} is undefined.
      
      The queue {Q} contains the sites that have been visited but
      are not yet finalized.  Apart from transients, every
      visited site is adjacent to a finalized site.  The {value}
      of each item {pv} in the queue is the difficulty of the easiest
      path found so far that leads from {pi} to {pv}, and is greater
      than or equal {rt[(*nt)-1]}.  The {value} of any of these nodes
      may still be revised downwards (but never below {rt[(*nt)-1]}).
    */
    
    /* Insert the starting node {pi} into the queue: */
    pqueue_insert(Q, pi, 0.0); (*nv)++;
    
    double rfmax = -INF; /* Maximum remoteness among the finalized nodes. */

    /* Dijkstra's loop: */
    while (pqueue_count(Q) > 0)
      {
        /* Pop the easiest site {pu} in queue, and its remoteness {ru}: */
        int pu = pqueue_head(Q); 
        double ru = pqueue_value(Q, pu);
        pqueue_delete(Q, pu);
        assert(ru >= rfmax); /* By Dijkstra's algorithm. */
        
        /* Are we done? */
        if (ru > dmax) { break; }
        
        /* Append it to the finalized sites: */
        site_list_append(L, pu, ru);

        /* Update {rfmax}: */
        rfmax = ru;
        
        /* Get {jini,jlim} so that {D->mat.e[jini..jlim-1]} are the edges out of {pu} in {D}: */
        int jini = D->start.e[pu]; 
        int jlim = D->start.e[pu+1];
        
        /* Visit the other endpoints of those edges: */
        int j;
        for (j = jini; j < jlim; j++)
          { /* Grab one edge {*ej} out of {pu}: */
            double_spmat_entry_t *ej = &(D->mat.e[j]);
            assert(ej->row == pu);
            
            /* Grab its destination {pv} and difficulty {dj}: */
            int pv = ej->col;
            double dj = ej->val;
            demand(dj >= 0, "step difficulty must be non-negative");
            
            /* Compute the total difficulty {rv} of the min path through {pu} to {pv}: */
            double rv = ru + dj;
            assert(rv >= rfmax);
            
            /* Check whether {pv} is finalized, visited but not finalized, or unvisited: */
            unsigned kv = site_list_position(L, pv);
            if (kv < L->n)
              { /* Site {pv} has already been finalized, ignore: */
                assert(rv >= L->val[kv]); /* Paranoia... */
              }
            else if (pqueue_has(Q, pv))
              { /* Site {pv} was visited but not finalized, update its value: */
                pqueue_set_value(Q, pv, rv);
              }
            else
              { /* Site {pv} is yet unvisited, insert it in {Q}: */
                pqueue_insert(Q, pv, rv); (*nv)++;
              }
          }
      }
      
    return;
  }

void compute_dir_distribution_params
  ( double xt, 
    double yt, 
    double dm, 
    double *xg, 
    double *yg, 
    double *dg
  );
  /* Computes the center {(xg,yg)} and the half-weight radius {dg}
    of the probability distribution of step directions,
    corresponding to a nonzero migration trend vector {(xt,yt)}
    and a Gaussian random drift component with half-weight
    radius {dm}.
    
    The procedure assumes that the desired probability distribution
    for the total trip displacement vector, scaled so that its length
    is the trip's difficulty {dt = hypot(xt,yt)}, is given by a
    Gaussian with half-weight radius {dm}, centered on {(xt,yt)} and
    restricted to the circle with radius {dt}. The distribution for
    step directions should therefore be this Gaussian, with the domain
    scaled down by {dt}.
    
    However, we need to consider that the step direction is quantized
    to 45 degree increments. Therefore, the radius {dg} is corrected
    to account by the corresonding quantization noise. */

distr_params_t *new_distr_params
  ( double shw_rel_diff, 
    double dpw_rel_diff, 
    double mean_trip_remt, 
    double start_site_wt,     
    bool_t use_new_state,      
    double mean_lang_diff,
    double vacant_factor
  )
  {
    distr_params_t *dpm = (distr_params_t *)notnull(malloc(sizeof(distr_params_t)), "out of mem"); 
    dpm->shw_rel_diff   = shw_rel_diff;
    dpm->dpw_rel_diff   = dpw_rel_diff;
    dpm->mean_trip_remt = mean_trip_remt;
    dpm->start_site_wt  = start_site_wt;
    dpm->use_new_state  = use_new_state;
    dpm->mean_lang_diff = mean_lang_diff;
    dpm->vacant_factor  = vacant_factor;  
    /* For now: */
    dpm->D = NULL;
    dpm->M = NULL;
    return dpm;
  }

void precompute_step_difficulties(geography_t *geo, distr_params_t *dpm)
  { 
    /* Check for redundant calls: */
    demand(dpm->D == NULL, "step diff matrix was already computed");
    
    /* Compute it: */
    dpm->D = compute_step_difficulty_matrix(geo, dpm->shw_rel_diff, dpm->dpw_rel_diff);
  }

void precompute_reachabilities(geography_t *geo, distr_params_t *dpm)
  {
    /* Check for redundant calls: */
    demand(dpm->M == NULL, "reachability matrix was already computed");
    
    /* We can precompute {dpm->M} only if remoteness is relevant: */
    demand(dpm->mean_trip_remt < +INF, "remoteness must be relevant");
    
    /* Make sure that the step difficulty matrix is present: */
    if (dpm->D == NULL) { precompute_step_difficulties(geo, dpm); }
    
    /* Choose the cutoff remoteness {dmax}. 
      We assume that, for any starting site {pi}, the weight of any 
      site {pk != pi} is bounded by {g(rk) = truncated_bell(rk/rm)},
      where {rk} is the remoteness of {pk} from {pi}, and
      {rm} is {dpm->mean_trip_remt}.  Therefore, the best cutoff for 
      {dpm->M} should be the smallest real {dmax} such that {g(dmax) == 0},
      namely {dmax == rm*TRUNCATED_BELL_MAX_Z}. So:
    */
    double dmax = dpm->mean_trip_remt * TRUNCATED_BELL_MAX_Z; 
    
    dpm->M = compute_reachability_matrix(geo, dpm->D, dmax);
  }

double compute_site_weight
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    site_id_t pk,        /* Id of site in question.  */ 
    double rk,           /* Its remoteness from the base site.  */ 
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t pi         /* Id of starting site.  */ 
  )
  {

    /* Reject sites on water: */
    if (get_altitude(geo, pk) <= 0) { return 0; }

    /* Compute weight {wk} as product of several factor: */
    double wk = 1;  /* For starters. */

    /* Apply the occupancy-related penalty factors: */
    if ((dpm->vacant_factor != 1.000) || (dpm->mean_lang_diff != INF))
      { /* Occupancy of site is significant: */
        frame_t *frm = (dpm->use_new_state ? new : old);
        site_t *ck = get_site_address(frm, pk);
        lang_t vk = ck->oc.lang;
        if (vk == NULL_LANG) 
          { /* Vacant site: */
            wk *= dpm->vacant_factor;
          }
        else if (dpm->mean_lang_diff != INF)
          { /* Gaussian function of language difference: */
            double zk = rel_lang_difference(ss->lang, vk)/dpm->mean_lang_diff; 
            double wlang = truncated_bell(zk);
            wk *= wlang;
          }
        if (wk == 0) { return 0; }
      }

    /* Apply self-avoidance penatly factors: */
    if (pk == pi)
      { /* Fixed weight factor: */
        wk = dpm->start_site_wt;
      }
    else
      { /* Apply the penalty fators related to remoteness: */
        if (dpm->mean_trip_remt == INF)
          { /* There is no penalty factor for remoteness: */ }
        else
          { /* Remoteness is relevant. */
            /* The penalty factor {wgeom} is a Gaussian function of {rk}: */
            double zk = rk/dpm->mean_trip_remt;
            double wgeom = truncated_bell(zk);
            wk *= wgeom;
          }
      }
    if (wk == 0) { return 0; }

    return wk;
  }

void collect_sites_from_mob_matrix
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t pi,        /* Site id of starting site.  */ 
    site_list_t *L,      /* (OUT) list of reachable sites. */
    int *mQ              /* (OUT) number of sites visited. */
  )
  {
    /* Get the reachability matrix {M}: */
    if (dpm->M == NULL) { precompute_reachabilities(geo, dpm); }
    mob_matrix_t *M = dpm->M;

    /* Get hold of edges out of {pi} in {M}: */
    int jini = M->start.e[pi];
    int jlim = M->start.e[pi+1];
    /* Edges are {M->mat.e[jini..jlim-1]}. */
    
    /* Flush the list {L}: */
    site_list_reset(L);

    /* Copy them to {L}, computing their probability weights: */
    int j;
    for (j = jini; j < jlim; j++)
      { /* Get the next reachable site {pj} and its remoteness {rj}: */
        double_spmat_entry_t *ej = &(M->mat.e[j]); /* Edge number {j}. */
        assert(ej->row == pi);
        site_id_t pj = ej->col; /* Its destination site. */
        double rj = ej->val;    /* The remoteness of that site. */
        /* Compute its probability weight: */
        double wj = compute_site_weight(old, new, geo, pj, rj, dpm, ss, pi);
        if (wj > 0)
          { /* Append {pj} to {Q}: */
            site_list_append(L, pj, rj); 
          }
      }

    /* Return the counts: */
    (*mQ) = jlim - jini;
 }

void collect_sites_dijkstra
  ( frame_t *old,        /* Old frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    double xt,           /* X coordinate of intended displacement. */
    double yt,           /* Y coordinate of intended displacement. */
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t pi,        /* Site id of starting site.  */ 
    site_list_t *L,      /* (OUT) list of reachable sites. */
    pqueue_t *Q,         /* (WORK) priority queue. */
    int *mQ              /* (OUT) number of sites visited. */
  )
  {
    int nx = geo->nx;
    /* int ny = geo->ny; */
    
    /* bool_t debug = ((*pp) == site_id(310,80, nx,ny)); */
    bool_t debug = FALSE;

    /* Get coordinates {(xi,yi)} of starting site: */
    int xi = site_X_from_id(pi,nx,ny);
    int yi = site_Y_from_id(pi,nx,ny);

    /* Select the frame {frm} that is relevant for weight computation: */
    frame_t *frm = (dpm->use_new_state ? new : old);
    
    auto double step_df(site_id_t pc, int xs, int ys);
       /* The difficulty of the step from {pc} to {pc+(xs,ys)}. */
    
    auto double remt_wt(double rc);
       /* An upper bound to the weight of a site with remoteness {rc}
         from {pi} (irrespective of its state). */
    
    auto double site_wt(double rc, site_id_t pc);
       /* Relative weight of site {pc} which has remoteness {rc}
         from {pi}. The actual weight is the result of this 
         function multiplied by {remt_wt(rc)}. */
    
    /* Gather sites by Dijkstra's algorithm: */
    double et = hypot(xt, yt); /* Length of trend vector. */
    double rt = et; /* Intended difficulty of trip in the direction of {(xt,yt)}. */
    generic_dijkstra(geo, pi, step_df, remt_wt, site_wt, L, Q, mQ);
    if (debug) 
      { print_site_list(stderr, frm, geo, "fixed sites\n", L, "\n");
        print_pqueue(stderr, frm, geo, "non-fixed sites\n", Q, "\n");
      }
    return;
    
    double step_df(site_id_t pc, int xs, int ys)
      { return step_difficulty(geo, pc, xs, ys, dpm->shw_rel_diff, dpm->dpw_rel_diff); }
    
    double remt_wt(double rc)
      { double wt;
        if (rc <= rt)
          { /* Must return 1 because of the non-decreasing requirement. */
            wt = 1;
          }
        else if (dpm->mean_trip_remt == INF)
          { /* There is no penalty for remoteness: */
            wt = 1;
          }
        else
          { /* Gaussian function of remoteness exceeding {rt}: */
            double ze = (rc-rt)/dpm->mean_trip_remt;
            wt = truncated_bell(ze);
          }
        return wt;
      }

    double site_wt(double rc, site_id_t pc)
      { 
        double wt = 1;  /* For starters. */

        /* Apply the occupancy-related penalty factors: */
        if ((dpm->vacant_factor != 1.000) || (dpm->mean_lang_diff != INF))
          { /* Occupancy of site is significant: */
            frame_t *frm = (dpm->use_new_state ? new : old);
            site_t *cc = get_site_address(frm, pc);
            lang_t vc = cc->oc.lang;
            if (vc == NULL_LANG) 
              { /* Vacant site: */
                wt *= dpm->vacant_factor;
              }
            else if (dpm->mean_lang_diff != INF)
              { /* Gaussian function of language difference: */
                double zc = rel_lang_difference(ss->lang, vc)/dpm->mean_lang_diff; 
                double wlang = truncated_bell(zc);
                wt *= wlang;
              }
            if (wt == 0) { return 0; }
          }
        
        /* Apply self-avoidance penatly factors: */
        if (pc == pi)
          { /* Fixed weight factor: */
            wt = dpm->start_site_wt;
          }
        else
          { /* Apply the penalty fators related to remoteness and direction: */
            if (rt == 0)
              { /* No trend vector -- the penalty factor was already supplied by {remt_wt}: */ }
            else if (dpm->mean_trip_remt == INF)
              { /* There is no penalty factor for remoteness: */ }
            else
              { /* The factor {remt_wt} was only an upper bound. */
                /* Get coordinates {(xc,yc)} of starting site: */
                int xc = site_X_from_id(pc,nx,ny);
                int yc = site_Y_from_id(pc,nx,ny);
                /* Compute the actual displacement {h = |(xh,yh)|} of the trip: */
                double xh = xc - xi;
                double yh = yc - yi;
                double h = hypot(xh, yh);
                demand(h > 0, "hypot bug");
                demand(rc > 0, "zero difficulty for non-zero displacement");
                /* Estimate the terrain difficulty factor {usc}: */
                double usc = rc/h;
                /* Scale the required displacement by the terrain difficulty: */
                double xu = xt/usc, yu = yt/usc;
                /* Estimate the remoteness {ru} between {(xu,yu)} and {(xh,yh)}: */
                double ru = hypot(xu - xh, yu - yh)*usc;
                /* The penalty factor {wgeom} is a Gaussian function of {ru}: */
                double wgeom;
                double zu = ru/dpm->mean_trip_remt;
                if (rc <= rt)
                  { /* {remt_wt} was 1, no correction needed: */
                    wgeom = truncated_bell(zu);
                  }
                else
                  { /* Compensate for {remt_wt}: */
                    double zm = (rc - rt)/dpm->mean_trip_remt;
                    wgeom = truncated_bell(zu)/truncated_bell(zm);
                  }
                wt *= wgeom;
                if (wt == 0) { return 0; }
              }
          }
          
        /* Make sure that weights are bounded by {remt_wt}: */
        if (wt > 1) { wt = 1; }
        return wt;
      }
  }

void select_site_from_mob_matrix
  ( frame_t *old,        /* Current frame. */
    frame_t *new,        /* New frame. */
    geography_t *geo,    /* World's geography. */
    distr_params_t *dpm, /* Parameters of probability distribution. */ 
    sibs_t *ss,          /* Attributes of selector. */ 
    site_id_t *pp,       /* (IN) Id of reference site; (OUT) Id of selected site. */ 
    site_list_t *L,      /* (OUT) list of reachable sites. */
    int *nQ,             /* (OUT) number of sites considered. */
    int *mQ              /* (OUT) number of sites visited. */
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    
    /* bool_t debug = ((*pp) == site_id(310,80, nx,ny)); */
    bool_t debug = FALSE;
    
    if (debug) 
      { fprintf(stderr, "\n");
        print_site(stderr, nx, ny, "    pick for", ss, pp, NULL, "\n");
      }
    
    /* Select the frame {frm} that is relevant for weight computation: */
    frame_t *frm = (dpm->use_new_state ? new : old);
    
    /* Save starting sie {pi}: */
    site_id_t pi = (*pp);
    
    /* Get the edges out of {pi} in {M}, and compute their probability weights: */
    collect_sites_from_mob_matrix(old, new, geo, dpm, ss, pi, L, mQ);
    (*nQ) = L->n;
     
    /* Select one site from the distribution: */
    int ks = pick_elem(L->val, L->n, 0, 0);
    site_id_t ps = L->itm[ks];

    if (debug) 
      { site_t *cs = &(frm->site[ps]);
        print_site(stderr, nx, ny, "    picked  ", &(cs->oc), &ps, NULL, "\n");
        fprintf(stderr, "\n");
      }

    /* Return its coordinates: */
    (*pp) = ps;
  }

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
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    
    /* bool_t debug = ((*pp) == site_id(310,80, nx,ny)); */
    bool_t debug = FALSE;
    
    if (debug) 
      { fprintf(stderr, "\n");
        print_site(stderr, nx, ny, "    pick for", ss, pp, NULL, "\n");
      }
    
    /* Select the frame {frm} that is relevant for weight computation: */
    frame_t *frm = (dpm->use_new_state ? new : old);
    
    /* Save starting sie {pi}: */
    site_id_t pi = (*pp);
    
    /* Collect relevant sites with this probability distribution: */
    collect_sites_dijkstra(old, new, geo, dpm, xt, yt, ss, pi, L, Q, mQ);
    (*nQ) = L->n;

    /* Select one site from the distribution: */
    int ks = pick_elem(L->val, L->n, 0, 0);
    site_id_t ps = L->itm[ks];

    if (debug) 
      { site_t *cs = &(frm->site[ps]);
        print_site(stderr, nx, ny, "    picked  ", &(cs->oc), &ps, NULL, "\n");
        fprintf(stderr, "\n");
      }

    /* Return its coordinates: */
    (*pp) = ps;

  }

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
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;

    bool_t debug = FALSE;
    
    /* Save the initial site {pi}: */
    site_id_t pi = (*pp);
    
    /* Select the frame {frm} that is relevant for weight computation: */
    frame_t *frm = (dpm->use_new_state ? new : old);
    
    /* This procedure is no good if the stop weight may be zero: */
    demand(dpm->vacant_factor > 0.05, "unsuitable {vacant_factor}");
    demand(dpm->mean_lang_diff > 0.01, "unsuitable {mean_lang_diff}");
    
    /* Compute the parameters {(xg,yg)} and {rg2} of the step direction distribution: */
    double xg, yg, dg;
    compute_dir_distribution_params(xt, yt, dpm->mean_trip_remt, &xg, &yg, &dg);
    
    /* Compute the desired trip difficulty {dt}: */
    double dt = hypot(xt,yt); 
    
    /* Perform a random walk governed by the given parameters. */
    
    auto double step_df(site_id_t pc, int xs, int ys);
      /* The difficulty of the step from {pc} to {pc + (xs,ys)}. */
    
    auto double halt_wt(double dc, site_id_t pc);
      /* A preference weight for halting the walk at the site {pc}
        which has remoteness {dc} from {pi}. */
    
    auto double step_wt(double es, double ds, site_id_t pc, int xs, int ys);
      /* A preference weight for taking the step from {pc} to
        {pc + (xs,ys)}, which has Euclidean length {es=hypot(xs,ys)}
        and difficulty {ds}. */
    
    random_walk(geo, pp, step_df, halt_wt, step_wt, L, nQ, mQ);
    if (debug) 
      { site_t *cp = &(frm->site[(*pp)]);
        print_site(stderr, nx, ny, "    picked  ", &(cp->oc), pp, NULL, "\n");
        fprintf(stderr, "\n");
      }
    return;
    
    double step_df(site_id_t pc, int xs, int ys)
      { return step_difficulty(geo, pc, xs, ys, dpm->shw_rel_diff, dpm->dpw_rel_diff); }
    
    double halt_wt(double dc, site_id_t pc)
      { 
        double wt = 1;  /* For starters. */
        
        /* Check whether site is on dry land: */
        if (get_altitude(geo, pc) <= 0) { return 0; }
        
        /* Apply the occupancy-related penalty factors: */
        if ((dpm->vacant_factor != 1.000) || (dpm->mean_lang_diff != INF))
          { /* Occupancy of site is significant: */
            frame_t *frm = (dpm->use_new_state ? new : old);
            site_t *cc = get_site_address(frm, pc);
            lang_t vc = cc->oc.lang;
            if (vc == NULL_LANG) 
              { /* Vacant site: */
                wt *= dpm->vacant_factor;
              }
            else if (dpm->mean_lang_diff != INF)
              { /* Gaussian function of language difference: */
                double zc = rel_lang_difference(ss->lang, vc)/dpm->mean_lang_diff; 
                double wlang = truncated_bell(zc);
                wt *= wlang;
              }
          }
        if (wt == 0) { return 0; }
        
        /* Apply self-avoidance penatly factors: */
        if (pc == pi)
          { /* Fixed weight factor: */
            wt = dpm->start_site_wt;
          }
        else
          { /* Apply the penalty fators related to position of final site: */
            if (dpm->mean_trip_remt == +INF)
              { /* There is no penalty factor for remoteness: */ }
            else
              { /* The geometric penalty factor depends only on the
                  total difficulty of the trip, since the direction of
                  the trip was (roughly) accounted for in
                  {step_wt}: */
                demand(dc > 0, "zero difficulty for non-zero displacement");
                /* Compute excess path difficulty {zc}, rel to {dpm->mean_trip_remt}: */
                double zc = (dc - dt)/dpm->mean_trip_remt;
                /* The penalty factor {wgeom} is a sigmoid function of {zc}: */
                double wgeom = truncated_sigmoid(zc/2);
                wt *= wgeom;
              }
          }
        if (wt == 0) { return 0; }
          
        /* Make sure that weights are bounded by {remt_wt}: */
        if (wt > 1) { wt = 1; }
        return wt;
      }
    
    double step_wt(double es, double ds, site_id_t pc, int xs, int ys)
      { 
        assert(es != 0); /* Staying put is not a valid step. */
        assert(ds != 0); /* A nonzero step must have a nonzero difficulty. */
        double wt;
        if (dpm->mean_trip_remt == INF)
          { /* There is no penalty for trip difficulty: */
            wt = 1.0;
          }
        else if (dt == 0)
          { /* No trend vector; weight decreases with difficulty: */
            double zs = ds/dpm->mean_trip_remt;
            wt = truncated_bell(zs);
          }
        else
          { /* Compute the unit vector {(xu,yu)} in the directon of the step: */
            double xu = xs/es;
            double yu = ys/es;
            /* Evaluate the Gaussian directional distribution on {(xg,yu)}: */
            double xh = xu - xg;
            double yh = yu - yg;
            double dh = hypot(xh, yh);
            double zh = dh/dg; 
            wt = truncated_bell(zh);
          }
        return wt;
      }
  }

void compute_dir_distribution_params
  ( double xt, 
    double yt, 
    double dm, 
    double *xg, 
    double *yg, 
    double *dg
  )
  { /* The center of the distrib is merely the direction of {(xt,yt)}: */
    double dt = hypot(xt,yt);
    (*xg) = xt/dt;
    (*yg) = yt/dt;
    /* The basic half-weight radius is {dm} scaled down by {dt}: */
    double db = dm/dt;
    /* The quantization error is uniform in {[-s _ +s]} where {s=sin(PI/8)}: */
    double dc = sin(M_PI/8)/sqrt(3);
    /* Combine the two: */
    (*dg) = hypot(db, dc); 
  }

void generic_dijkstra
  ( geography_t *geo,           /* World's geography. */
    site_id_t pi,               /* Site id of starting site. */
    step_diff_func_t step_df,   /* Defines the difficulty of a step. */
    remt_weight_funt_t remt_wt, /* Defines the remoteness-related weight of a site. */
    site_weight_func_t site_wt, /* Defines the endsite-related weight factor. */
    site_list_t *L,             /* (OUT) list of reachable sites. */
    pqueue_t *Q,                /* (WORK) priority queue. */
    int *mQ                     /* (OUT) number of sites visited. */
  )
  { 
    int nx = geo->nx;
    int ny = geo->ny;
    bool_t debug = FALSE;
    
    /* Flush {L} and {Q}: */
    site_list_reset(L);
    pqueue_reset(Q);
    pqueue_set_order(Q, +1);

    /* Sites in the list {L} and in the priority queue {Q} are the
      /visited/ sites, namely those for which we have already found a
      path from {pi} with finite difficulty. The value associated to
      each site in {Q} or {L} is the difficulty of this path.
      
      The lst {L} contains the visited sites which have been
      /finalized/, meaning that their value (remoteness from {pi}) is
      definitive. The list is sorted by increasing value.
      
      The queue {Q} is arranged as a heap, so that the element with
      lowest value is at position 0. The other entries are not yet
      finalized; their value may decrease and their positions may
      change. However, no value in {Q} can become smaller than any
      value in {L}. */
    
    int nv = 0; /* Count sites seen. */
    
    auto void enqueue_and_mark(site_id_t pf, double df);
      /* If {pf} is in {L}, the operation is a no-op.
        
        If {pf} is in {Q} and {df} is greater than or equal to its
        current value, the operation is a no-op.
        
        If {pf} is in {Q} and {df} is strictly less than the current
        value, updates its value to {df}, and rearranges {Q}
        accordingly.
        
        If {pf} is neither in {Q} nor in {L}, inserts it in {Q} with
        value {df}. */
    
    /* Initialize queue with starting site {(x,y)}: */
    enqueue_and_mark(pi, 0.0);
    
    double wtot = 0.0; /* Total weight of finalized sites. */

    /* Dijkstra's loop: */
    while (pqueue_count(Q) > 0)
      {
        /* Get the easiest site {pq} from {Q}, and its remoteness {dq}: */
        site_id_t pq = pqueue_head(Q);
        double dq = pqueue_value(Q, pq);
        pqueue_delete(Q, pq);
        /* Any still-unseen site has remoteness {dq} or more, so {pq} is final. */
        
        /* Compute weight factor {wremt} due to trip remoteness: */
        double wremt = remt_wt(dq);
        demand(wremt >= 0, "bad remt_wt()");
        
        /* Compute weight factor {wsite} due to other factors: */
        double wsite = site_wt(dq, pq);
        demand((wsite >= 0) && (wsite <= 1), "bad site_wt()");
          
        /* Compute the full weight: */
        double wq = wremt*wsite;
          
        /* Make {pq} final, with value {wq} (instead of {dq}): */
        site_list_append(L, pq, wq);
        
        if (debug)
          { print_site_coords(stderr, nx, ny, "      fixed site ", pq, ""); 
            fprintf(stderr, " remt = %8.3f", dq); 
            fprintf(stderr, " wremt = %6.4f wsite = %6.4f w = %6.4f\n", wremt, wsite, wq);
          }
        
        /* Update the total weight {wtot} of the fixed sites: */
        wtot += wq;
        
        /* Estimate the max total weight of unfixed sites: */
        double wrem = (nx*ny - nv)*wremt;
        
        /* If we got most of the available weight, stop: */
        if (wrem < 0.001*wtot) { break; }
        
        /* Check sites adjacent to {pq}: */
        int xs, ys;
        for (xs = -1; xs <= +1; xs++)
          { for (ys = -1; ys <= +1; ys++)
              { /* Get the neighbor's Y coord {yc}: */
                if ((xs != 0) || (ys != 0))
                  { /* Compute the difficulty of step {(xs,ys)} from {pq}: */
                    double ds = step_df(pq, xs, ys);
                    if (ds < +INF)
                      { /* Compute the coordinates {(xc,yc)} of the step's destination: */
                        int xc = site_X_from_id(pq,nx,ny) + xs;
                        while (xc < 0) { xc += nx; }
                        while (xc >= nx) { xc -= nx; }
                        int yc = site_Y_from_id(pq,nx,ny) + ys;
                        assert((yc >= 0) && (yc < ny)); 
                        /* Get the id {pc} of the step's destination: */
                        site_id_t pc = site_id(xc, yc, nx, ny);
                        /* Compute the difficulty of this path from start to {pc}: */
                        double dc = dq + ds;
                        /* Add/update remoteness of {pc} and re-sort the queue: */
                        enqueue_and_mark(pc, dc);
                      }
                  }
              }
          }
      }
    
    /* Return the number of visited sites in {*mQ}: */
    (*mQ) = nv;
    return;
    
    void enqueue_and_mark(site_id_t pf, double df)
      { assert(df < +INF); /* Otherwise we shouldn't have called it. */
        if (site_list_has(L, pf))
          { /* The site {pf} was alrady fixed: */ }
        else if (pqueue_has(Q, pf))
          { /* We have seen this site before: */
            if (df < pqueue_value(Q, pf))
              { /* This path is better than the previously found one; */
                pqueue_set_value(Q, pf, df);
                if (debug)
                  { print_site_coords(stderr, nx, ny, "        reached again ", pf, "");
                    fprintf(stderr, " with difficulty = %8.3f\n", df);
                  }
              }
          }
        else
          { /* First visit to this site: */
            pqueue_insert(Q, pf, df);
            if (debug)
              { print_site_coords(stderr, nx, ny, "        reached first ", pf, "");
                fprintf(stderr, " with difficulty = %8.3f\n", df);
              }
            nv++;
          }
      }
  }

void random_walk
  ( geography_t *geo,           /* World's geography. */
    site_id_t *pp,              /* (IN/OUT) Site id of initial and final site. */
    step_diff_func_t step_df,   /* Defines the difficulty of a step. */
    halt_weight_func_t halt_wt, /* Defines the probability weight of halting at a site. */
    step_weight_funt_t step_wt, /* Defines the probability weight of taking a step. */
    site_list_t *L,             /* Work area. */
    int *nQ,                    /* (OUT) Number of steps taken along the walk. */
    int *mQ                     /* (OUT) number of sites looked at along the walk. */
  )
  { 
    int nx = geo->nx;
    int ny = geo->ny;
    
    bool_t debug = ((*pp) == site_id(310,80, nx,ny));
    /* bool_t debug = FALSE; */
    
    int n_sites = 0; /* Number of sites in walk. */
    int n_steps = 0; /* Number os steps considered along the walk. */
    
    /* Grab the starting site {pi}: */
    site_id_t pi = (*pp);
    demand(pi < nx*ny, "bad {pi}");
    
    /* Current end of the walk: */
    site_id_t pc = pi;
    
    if (debug) { print_site_coords(stderr, nx, ny, "    walk from ", pc, "\n"); }

    /* Current difficulty of the trip: */
    double dc = 0;
    
    /* Number of alternatives at each site: */
    int na = 1 + 8; /* Halt here, or take a king step. */
    
    /* Probability weights, displacements, and difficulties of each alternative: */
    double wa[na], xa[na], ya[na], da[na];
    
    /* Main loop: */
    while (TRUE)
      { 
        /* One more site in the walk: */
        n_sites++; 
        
        if (debug)
          { print_site_coords(stderr, nx, ny, "      site at ", pc, "");
            fprintf(stderr, " diff = %8.3f\n", dc);
          }
        
        /* Get the weights {w[0..ns-1]} for valid steps: */
        int ia = 0; /* Index of next alternative. */
        int xs, ys;
        for (xs = -1; xs <= +1; xs++)
          { for (ys = -1; ys <= +1; ys++)
              { /* One more step considered: */
                n_steps++; 
                /* Compute the difficulty {ds} of this alternative: */
                double ds; 
                if ((xs == 0) && (ys == 0))
                  { /* Halting here is no sweat: */ ds = 0; }
                else
                  { /* Compute the difficulty of step {(xs,ys)} from {pc}: */
                    ds = step_df(pc, xs, ys);
                  }
                if (debug)
                  { fprintf(stderr, "        step (%+02d %+02d)  diff = %8.3f", xs, ys, ds); }

                /* Compute the weight {ws} of this alternative: */
                double ws; 
                if ((xs == 0) && (ys == 0))
                  { /* The alternative is to halt here: */
                    ws = halt_wt(dc, pc); 
                  }
                else
                  { /* Compute the probability weight of that step: */
                    if (ds >= +INF)
                      { /* Step is impossible, do not take: */ ws = 0; }
                    else
                      { /* Compute the Euclidean length of this step: */
                        double es = abs(xs) + abs(ys);
                        if (es == 2.0) { es = M_SQRT2; }
                        /* Compute the probability weight of this step: */
                        ws = step_wt(es, ds, pc, xs, ys);
                      }
                  }
                if (debug) { fprintf(stderr, "  whgt = %8.3f\n", ws); }

                /* Append {ws,xs,ys,ds} to the alternative table: */
                wa[ia] = ws; xa[ia] = xs; ya[ia] = ys; da[ia] = ds;
                ia++;
              }
          }
        assert(ia == na); 
        
        /* Choose one of the alternatives: */
        int ka = pick_elem(wa, na, -1, 0);
        /* Set {(xk,yk)} to its step, {dk} to its difficulty: */
        int xk, yk; 
        double dk;
        if (ka < 0)
          { /* Pick failed, assume halt: */ xk = 0; yk = 0; dk = 0; }
        else
          { /* Pick succeeded: */ xk = xa[ka]; yk = ya[ka]; dk = da[ka]; }
          
        if (debug)
          { fprintf(stderr, "      chosen (%+02d %+02d)  diff = %8.3f\n", xk, yk, dk); }

        /* Have we decided to stop? */
        if ((xk == 0) && (yk == 0)) { /* Exit the main loop: */ break; }
          
        /* Compute the coords {(xd,yd)} of the destination site: */
        int xd = site_X_from_id(pc,nx,ny) + xk;
        while (xd < 0) { xd += nx; }
        while (xd >= nx) { xd -= nx; }
        int yd = site_Y_from_id(pc,nx,ny) + yk;
        assert((yd >= 0) && (yd < ny));
        /* Move {pc} the destination site: */
        pc = site_id(xd,yd, nx,ny);
        
        /* Update the total difficulty of the walk so far: */
        dc += dk;
      }

    /* Return the final site: */
    (*pp) = pc; 
    /* Return the site counts: */
    (*nQ) = n_sites;
    (*mQ) = n_steps;
    return;
  }
  
void print_site_list
  ( FILE *wr, 
    frame_t *frm, 
    geography_t *geo, 
    char *hd, 
    site_list_t *L,
    char *ft 
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    
    if ((hd != NULL) && ((*hd) != 0)) { fputs(hd, wr); } 
    int k;
    for (k = 0; k < L->n; k++)
      { site_id_t pk = L->itm[k];
        double wk = L->val[k]; /* Weight. */
        site_t *ck = &(frm->site[pk]);
        double alt = get_altitude(geo, pk);
        fprintf(wr, "      cand %4d   wgth = %7.4f", k, wk);
        print_site(wr, nx, ny, " ", &(ck->oc), &pk, &alt, "\n");
      }
    if ((ft != NULL) && ((*ft) != 0)) { fputs(ft, wr); } 
  }
 
  
void print_pqueue
  ( FILE *wr, 
    frame_t *frm, 
    geography_t *geo, 
    char *hd, 
    pqueue_t *Q, 
    char *ft 
  )
  {
    int nx = geo->nx;
    int ny = geo->ny;
    
    if ((hd != NULL) && ((*hd) != 0)) { fputs(hd, wr); } 
    int k;
    int n = pqueue_count(Q);
    for (k = 0; k < n; k++)
      { site_id_t pk = pqueue_item(Q, k);
        double rk = pqueue_value(Q, pk); /* Remoteness. */
        site_t *ck = &(frm->site[pk]);
        double alt = get_altitude(geo, pk);
        fprintf(wr, "      cand %4d   remt = %7.4f", k, rk);
        print_site(wr, nx, ny, " ", &(ck->oc), &pk, &alt, "\n");
      }
    if ((ft != NULL) && ((*ft) != 0)) { fputs(ft, wr); } 
  }
