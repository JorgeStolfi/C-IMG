/* Site-value tables with fast lookup. */
/* Last edited on 2024-12-21 14:00:19 by stolfi */ 

#ifndef langev_list_H
#define langev_list_H

#include <stdio.h>

#include <bool.h>

#include <langev_geog.h>

/* !!! Deserves to become a generic library routine. */

typedef struct site_list_t
  { int nsites;       /* Count of sites in the world. */
    int n;            /* Count of sites in list. */
    site_id_t *itm;   /* Items (site ids). */
    double *val;      /* Values associated to items. */
    unsigned *pos;    /* Maps items to positions in {itm}. */
  } site_list_t;
  /* A {site_list_t} is a list of pairs {(p,v)} where {p} is a site id
    and {v} is a real value associated to {p}.
    
    The field {nsites} is the number of sites in the world, so the
    site ids are numbers in the range {0..nsites-1}. The three vectors
    {itm,val,pos} have {nsites} elements.
    
    At any time, the ids of the sites in the list are {itm[0..n-1]},
    and they must be all distinct. The value associated to site
    {itm[k]}, for each {k} in {0..n-1 is}, is {val[k]}. A site with id
    {p} is stored in the list if and only if {pos[p]} is in {0..n-1},
    and {itm[pos[p]] == p}.
    
    The order of the pairs in the list is modified only by the 
    procedures {site_list_swap} and {site_list_delete}. 
    
    The client must not modify any field of {L} except through the
    procedures in this module.  All times below are worst-case. */

site_list_t *site_list_new(int nsites);
  /* Returns a pointer to a newly allocated {site_list_t} record,
    with internal tables allocated for a world with {nsites} sites. */

unsigned site_list_count(site_list_t *L);
  /* The number of items currently in {L}. Time: {O(1)}.  */

bool_t site_list_has(site_list_t *L, site_id_t p);
  /* TRUE iff item {p} is in {L}. Time: {O(1)}. */

unsigned site_list_position(site_list_t *L, site_id_t p);
  /* Returns the index {ix} in {0..Q.n-1} such that {Q->el[ix] == z};
     or {Q.n} if no such {ix} exists. Time: {O(1)}. */

double site_list_value(site_list_t *L, site_id_t p);
  /* Current value of item {p} in {L}; undefined if {p} is not
    in list. Time: {O(1)}.  */

double site_list_item(site_list_t *L, unsigned i);
  /* The item currently at position {i} of the list, where {i}
    may range from 0 to {site_list_count(L)-1}. Time: {O(1)}. */

/* MODIFYING */

void site_list_append(site_list_t *L, site_id_t p, double v);
  /* Appends the item {p} to {L}, with value {v}. Undefined if {p} is
    already in {L}. Does not change the positions of pre-existing
    elements. Time: {O(1)} worst-case if no storage re-allocation is
    needed; otherwise {O(n)} worst-case, {O(1)} amortized. */

void site_list_delete(site_list_t *L, site_id_t p);
  /* Removes the item {p} from {L}. Undefined if {p}
    is not in {L}.  If {p} is not the last element,
    alsomoves the last element to the position of {p}.
    Time: {O(1)}. */

void site_list_swap(site_list_t *L, unsigned i, unsigned j);
  /* Swaps the entries at positions {i} and {j} of {L}. 
    The two indices must be in the range {0..Q.n-1}. Time: {O(1)}. */

void site_list_set_value(site_list_t *L, site_id_t p, double v);
  /* Sets the value of item {p} to {v}. Time: {O(1)}. */

void site_list_reset(site_list_t *L);
  /* Deletes all elements.  Time: {O(1)}. Does not free any internal 
    storage. */

void site_list_free(site_list_t *L);
  /* Frees the internal storage, and {*L} itself. */

int site_list_dblcmp(double x, double y);
  /* Returns -1, 0, or +1 depending on whether {x} is less than,
    equal to, or greater than {y}. */

void site_list_check(site_list_t *L);
   /* Checks whether the invariants of {L} are satisfied. */
   
#endif
