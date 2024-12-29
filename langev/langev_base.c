/* See {langev_move.h} */
/* Last edited on 2024-12-21 14:00:40 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <values.h>

#include <affirm.h>
#include <bool.h>
#include <vec.h>
#include <jspnm.h>
#include <jsfile.h>
#include <uint16_image.h>
#include <rn.h>

#include <langev_lang.h>
#include <langev_gene.h>
#include <langev_name.h>
#include <langev_base.h>
#include <langev_geog.h>
#include <langev_util.h>

frame_t *new_frame(int nx, int ny)
  { frame_t *frm = (frame_t*)notnull(malloc(sizeof(frame_t)), "out of memory");
    frm->cols = nx; frm->rows = ny;
    /* Allocate site array: */
    frm->site = (site_t *)notnull(malloc(ny*nx*sizeof(site_t)), "out of memory");
    /* Allocate language-related attribute tables: */
    frm->prs = rn_alloc(N_GR_LANG);
    frm->xtv = rn_alloc(N_GR_LANG);
    frm->ytv = rn_alloc(N_GR_LANG);
    return frm;
  }    

site_t *get_site_address(frame_t *frm, site_id_t p)
  { demand(p < frm->cols*frm->rows, "invalid site id");
    return &(frm->site[p]); 
  }

bool_t do_settle_twins
  ( frame_t *frm,
    geography_t *geo,   /* World's geography. */
    sibs_t *ss,
    site_id_t pp, 
    bool_t verbose
  )
  { 
    int nx = geo->nx;
    int ny = geo->ny;

    double altitude = get_altitude(geo, pp);
    if (verbose) { print_site(stderr, nx, ny, "settling twins ", ss, &pp, &altitude, "\n"); }
    if (altitude <= 0)
      { /* Site is on water -- the twins die: */
        if (verbose) { fprintf(stderr, "twins drowned\n"); }
        return FALSE;
      }
    else
      { /* Site is on dry land -- get it: */
        site_t *c = get_site_address(frm, pp);
        if (verbose) { print_site(stderr, nx, ny, "present status ", &(c->oc), NULL, NULL, "\n"); }
        /* Kill any previous occupants and put {*s} in its place: */
        c->oc = (*ss);
        return TRUE;
      }
  }

void fill_sites(frame_t *frm, sibs_t sibs, double remt)
  { int nx = frm->cols;
    int ny = frm->rows;
    int nsites = nx*ny;
    site_id_t p;
    for (p = 0; p < nsites; p++)
      { site_t *c = &(frm->site[p]);
        c->oc = sibs; c->remt = remt;
      }
  }

double get_language_prestige(frame_t *frm, lang_t v)
  { /* The current implementation simply uses the language's major group: */
    int k = language_group(v);
    assert((k >= 0) && (k < N_GR_LANG));
    return frm->prs[k];
  }

void get_language_drift_vector(frame_t *frm, lang_t v, double *xt, double *yt)
  { /* The current implementation simply uses the language's major group: */
    int k = language_group(v);
    assert((k >= 0) && (k < N_GR_LANG));
    (*xt) = frm->xtv[k];
    (*yt) = frm->ytv[k];
  }

void output_frame
  ( char *prefix, 
    int iep, 
    frame_t *frm, 
    geography_t *geo
  )
  { 
    int nx = frm->cols;
    int ny = frm->rows;
    
    bool_t debug = FALSE;
    
    auto void lang_color(site_id_t pc, sibs_t *sc, float rgb[]);
      /* Stores into {smp[0..nc-1]} an RGB color that depends on the
        language component of {*sc}, assumed to be the atrtibutes of
        the current occupants of the frame {frm}. If the site is
        vacant, however, the color will depend on its altitude. */
  
    void lang_color(site_id_t pc, sibs_t *sc, float rgb[])   
      { /* Get relative altitude in range {[-1 _ +1]}: */
        double altitude = get_altitude(geo, pc);
        /* Paranoia: */
        demand((altitude > 0) || (sc->lang == NULL_LANG), "site on water is not vacant");
        /* map_altitude_to_color(altitude, rgb); */
        if (sc->lang == NULL_LANG)
          { /* Vacant site -- display the altitude as black/grey/white: */
            double Y;
            if (altitude < 0)
              { /* Deep water: */ Y = 0.600; }
            else if (altitude > 0)
              { /* Dry land: */ Y = 0.950; }
            else
              { /* Shallow water: */ Y = 0.800; }
            rgb[0] = rgb[1] = rgb[2] = Y;
          }
        else
          { 
            if (debug) { print_site(stderr, nx, ny, "alive ", sc, &pc, &altitude, "\n"); }
            /* Occupied dite - display language color: */
            map_language_to_color(sc->lang, rgb);
          }
      }

    uint16_image_t *img = colorize_frame(frm, lang_color);
    char *fileName = jsprintf("%s-%06d.ppm", prefix, iep);
    FILE *wr = open_write(fileName, TRUE);
    bool_t forceplain = FALSE;
    bool_t verbose = TRUE;
    uint16_image_write_pnm_file(wr, img, forceplain, verbose);
    fclose(wr);
    uint16_image_free(img);
    free(fileName);
  } 

void print_site(FILE *wr, int nx, int ny, char *pre, sibs_t *sp, site_id_t *pp, double *altp, char *suf)
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    if (sp != NULL) { print_sibs(wr, " ", sp, ""); }
    if (pp != NULL) { print_site_coords(wr, nx, ny, " at ", *pp, ""); }
    if (altp != NULL) { print_altitude(wr, " alt = ", *altp, ""); }
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void print_coord(FILE *wr, char *pre, int z, char *suf)
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%4d", z);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void print_sibs(FILE *wr, char *pre, sibs_t *sp, char *suf)
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    print_name(wr, "", sp->name, "");
    print_name(wr, " = (", sp->ndad, "");
    print_name(wr, " + ", sp->nmom, ")");
    print_lang(wr, " lang = ", sp->lang, "");
    print_gene(wr, " gene = ", sp->gene, "");
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void print_site_coords(FILE *wr, int nx, int ny, char *pre, site_id_t p, char *suf)
  { 
    int x = site_X_from_id(p, nx, ny);
    int y = site_Y_from_id(p, nx, ny);
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    print_coords(wr, "(", x, " ", y, ")");
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void print_site_id(FILE *wr, char *pre, site_id_t p, char *suf)
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%8d", p);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void print_coords(FILE *wr, char *pre, int x, char *sep, int y, char *suf)
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%4d", x);
    if ((sep != NULL) && ((*sep) != 0)) { fputs(sep, wr); }
    fprintf(wr, "%4d", y);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void print_altitude(FILE *wr, char *pre, double alt, char *suf)
  { 
    if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fprintf(wr, "%6.3f", alt);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }
