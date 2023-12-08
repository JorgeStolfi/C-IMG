/* See langev_list.h */
/* Last edited on 2023-02-25 16:04:37 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <rn.h>

#include <langev_list.h>

site_list_t *site_list_new(int nsites)
  { site_list_t *L = (site_list_t *)notnull(malloc(sizeof(site_list_t)), "out of memory");
    L->nsites = nsites;
    L->n = 0;
    L->itm = (unsigned *)notnull(malloc(nsites*sizeof(unsigned)), "no mem"); 
    L->val = rn_alloc(nsites); 
    L->pos = (unsigned *)notnull(malloc(nsites*sizeof(unsigned)), "no mem");
    return L;
  }

void site_list_free(site_list_t *L)
  { free(L->itm);
    free(L->val);
    free(L->pos);
    free(L);
  }

/* QUERIES (WITHOUT SIDE EFFECTS) */

unsigned site_list_count(site_list_t *L)
  { return L->n; }

unsigned site_list_position(site_list_t *L, site_id_t z)
  { if (z >= L->nsites) { return L->n; }
    unsigned i = L->pos[z];
    if ((i > L->n) || (L->itm[i] != z)) { return L->n; }
    return i;
  }

bool_t site_list_has(site_list_t *L, site_id_t z)
  { unsigned i = site_list_position(L, z);
    return (i < L->n);
  }

double site_list_value(site_list_t *L, site_id_t z)
  { unsigned i = site_list_position(L, z);
    demand(i < L->n, "item not in queue");
    return L->val[i];
  }

double site_list_item(site_list_t *L, unsigned i)
  { demand(i < L->n, "no such position in queue");
    return L->itm[i];
  }

/* MODIFYING */

void site_list_append(site_list_t *L, site_id_t z, double v)
  { unsigned i = site_list_position(L, z);
    demand(i == L->n, "item already in queue");
    L->itm[i] = z;
    L->val[i] = v;
    L->pos[z] = i;
    L->n++;
  }

void site_list_delete(site_list_t *L, site_id_t z)
  { unsigned i = site_list_position(L, z);
    demand(i < L->n, "item not in queue");
    /* Make sure that element {z} is the last one: */
    unsigned j = L->n - 1; /* Index of last slot. */
    if (i != j) { site_list_swap(L, i, j); }
    /* Delete {z} assuming it is the last element: */
    L->n --; 
  }

void site_list_swap(site_list_t *L, unsigned i, unsigned j)
  {
    demand(i < L->n, "invalid index {i}");
    demand(j < L->n, "invalid index {j}");
    if (i != j)
      { /* Swap values and items: */
        site_id_t pi = L->itm[i];
        site_id_t pj = L->itm[j];
        L->itm[i] = pj;
        L->itm[j] = pi; 
        double vi = L->val[i]; 
        double vj = L->val[j]; 
        L->val[i] = vj;
        L->val[j] = vi;
        /* Update positions: */
        L->pos[pi] = j;
        L->pos[pj] = i;
      }
  }

void site_list_set_value(site_list_t *L, site_id_t z, double v)
  { unsigned i = site_list_position(L, z);
    demand(i < L->n, "item not in queue");
    L->val[i] = v;
  }

void site_list_reset(site_list_t *L)
  { L->n = 0; }

void site_list_check(site_list_t *L)
  {
    /* Validate {L->n}: */
    assert((L->n >= 0) && (L->n <= L->nsites)); 
    unsigned i;
    for (i = 0; i < L->n; i++)
      { /* Get item {z} in slot {i}: */
        site_id_t z = L->itm[i];
        /* Items must be in range {0..L->nsites-1}: */
        assert(z < L->nsites); 
        /* Get claimed position of item {z}: */
        unsigned k = L->pos[z];
        /* Must be this position: */
        assert(k == i);
      }
  }
  
int site_list_dblcmp(double x, double y) 
  { return (x < y ? -1 : ( x > y ? +1 : 0)); }
