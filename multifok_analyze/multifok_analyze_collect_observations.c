/* See {multifok_analyze_collect_observations.h}. */
/* Last edited on 2024-12-21 13:59:40 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsstring.h>
#include <rn.h>
#include <float_image.h>
#include <multifok_focus_op.h>

#include <multifok_analyze_extract_window.h>
#include <multifok_analyze_compute_rel_similarities.h>

#include <multifok_analyze_collect_observations.h>

void multifok_analyze_collect_observations
  ( float_image_t *ref,
    int32_t NF,
    float_image_t **frame, 
    int32_t NW, 
    double noise,
    bool_t verbose,
    int32_t *NDP,
    int32_t **ixP,
    int32_t **iyP,
    int32_t *NPP,
    double **e2P,
    int32_t *NQP, 
    char ***tnameP,
    double **termsP, 
    float_image_t ***qtimgP
  )
  {
    int32_t NC = (int32_t)ref->sz[0];
    int32_t NX = (int32_t)ref->sz[1];
    int32_t NY = (int32_t)ref->sz[2];
    demand(NC == 1, "images must be monochromatic");
    
    demand((NW >= 3) && ((NW % 2) == 1), "invalid window size");
    int32_t HW = (NW-1)/2; /* Window radius. */

    /* Compute the window basis: */
    if (verbose) { fprintf(stderr, "computing focus operator tables...\n"); }
    int32_t NS = NW*NW; /* Number of pixels (samples) in each window. */
    double *w = multifok_focus_op_prod_weights(NW); /* Length {NS}. */
    double **phi = multifok_focus_op_basis(NW); /* {NS} vectors of {NS} elements. */
    if (verbose) { multifok_focus_op_basis_check(NW, phi, w); }
    
    /* The window samples are normalized and then transformed to an orthonormal basis {phi}. */
    /* The quadratic terms are constructed from these transformed samples. */
    
    auto char *describe_term_3x3(int32_t q);
      /* Returns a string that describes the quadratic term with index {q},
        for {q} in {0..NQ-1}, for a {3x3} window. 
        If {q} is not in {0..NQ-1}, returns {NULL}.  */
        
    auto double compute_term_3x3(int32_t ns, double vr[], int32_t q);
      /* Returns the quadratic term with index {q}, for {q} in {0..NQ-1},
      given the normalized and remapped window contents {vr[0..ns-1]}.
      Specific for a {3x3} window. Blows up if {q} is not in {0..NQ-1}. */
      
    /* Discover the number of quadratic terms and their names: */
    bool_t cross_terms = FALSE; /* Include cross terms {P[i]*P[j]} in the observations, else only quadratics. */
    int32_t NQmax = NS*(NS+1)/2; /* At most each product of window ops is a term. */
    demand(NW == 3, "only 3x3 window is implemented for now");
    char **tname = notnull(malloc(NQmax*sizeof(char*)), "no mem"); /* To be trimmed later. */
    int32_t NQ = 0;
    while(TRUE)
      { char *tnq = describe_term_3x3(NQ);
        if (tnq == NULL) { break; }
        assert(NQ < NQmax);
        tname[NQ] = tnq;
        NQ++;
      }
    if (verbose) { fprintf(stderr, "generating %d quadratic terms per observation\n", NQ); }
    tname = notnull(realloc(tname, NQ*sizeof(char*)), "no mem");  
    
    /* Allocate the quadratic term images, fill them with zeros: */
    float_image_t **qtimg = notnull(malloc(NF*NQ*sizeof(float_image_t *)), "no mem");
    for (uint32_t f = 0;  f < NF; f++)
      { for (uint32_t q = 0;  q < NQ; q++)
          { float_image_t *img = float_image_new(1, NX, NY); 
            float_image_fill(img, 0.0f);
            qtimg[f*NQ + q] = img;
          }
      }
  
    /* Collect the observations: */
    int32_t NDmax = (NX - NW + 1)*(NY - NW + 1); /* Max number of window positions in each image. */
    int32_t NPmax = NDmax*NF; /* Max number of observations extracted. */
    int32_t *ix = notnull(malloc(NDmax*sizeof(int32_t)), "no mem");
    int32_t *iy = notnull(malloc(NDmax*sizeof(int32_t)), "no mem");
    double *e2 = notnull(malloc(NPmax*sizeof(double)), "no mem");
    double *terms = notnull(malloc(NPmax*NQ*sizeof(double)), "no mem");
    double dev_min = 3*noise; /* Min sample deviation to consider window significant. */
    
    double rr[NS]; /* Window contents from ref image, normalized. */
    double fr[NS]; /* Window contents from frame image, normalized. */
    double gr[NS]; /* Orthogonal transform of {fr[0..NS-1]}. */
    
    int32_t ND = 0; /* Number of significant window positions in ref image found so far. */
    int32_t NP = 0; /* Number of observations collected so far. */
    int32_t nstrange = 0; /* Number of frame windows with more energy than the reference window. */
    for (int32_t dx = HW; dx < NX - HW; dx++)
      { for (int32_t dy = HW; dy < NY - HW; dy++)
          { /* Extract samples {rr[0..NW-1]} from ref image in window centered at col {dx} and row {dy}: */
            double rdev = multifok_analyze_extract_window(NW, ref, dx, dy, noise, rr, w);
            if (rdev >= dev_min)
              { assert(ND < NDmax);
                int32_t d = ND; /* Index of this relevant window position. */
                /* Store the window position: */
                ix[d] = dx; iy[d] = dy;
                /* Extract observations from frame images: */
                for (uint32_t f = 0;  f < NF; f++)
                  { assert(NP < NPmax);
                    int32_t p = NP;  /* Index of this observation. */
                    assert(p == d*NF + f);
                    /* Extract samples {fr[0..NW-1]} from frame {f} in window centered ar {dx,dy}: */
                    float_image_t *fimg = frame[f];
                    double fdev = multifok_analyze_extract_window(NW, fimg, dx, dy, noise, fr, w);
                    if (fdev > 1.5*rdev) { nstrange++; }
                    /* Compute the dependent variable: */
                    double D2 = multifok_focus_op_dist_sqr(NW, rr, fr, w);
                    e2[p] = 1.0 - D2/4;
                    /* Transform {fr} to {gr}: */
                    multifok_focus_op_remap_samples(NW, fr, phi, w, gr);
                    double *ts = &(terms[p*NQ]);  /* Quadratic terms from this frame window. */
                    /* Generate the quadratic terms from remapped samples, save in images: */
                    for (uint32_t q = 0;  q < NQ; q++) 
                      { ts[q] = compute_term_3x3(NS, gr, q);
                        float_image_t *img = qtimg[f*NQ + q];
                        float_image_set_sample(img, 0, dx, dy, (float)(ts[q]));
                      }
                    NP++;
                  }
                ND++;
              }
          }
      }
    if (nstrange > 0) { fprintf(stderr, "found %d anomalously energetic frame windows\n", nstrange); }
    /* Trim observation arrays: */
    ix = notnull(realloc(ix, ND*sizeof(int32_t)), "no mem");
    iy = notnull(realloc(iy, ND*sizeof(int32_t)), "no mem");
    e2 = notnull(realloc(e2, NP*sizeof(double)), "no mem");
    terms = notnull(realloc(terms, NP*NQ*sizeof(double)), "no mem");

    (*NDP) = ND;
    (*ixP) = ix;
    (*iyP) = iy;
    (*NPP) = NP;
    (*e2P) = e2;
    (*NQP) = NQ;
    (*tnameP) = tname;
    (*termsP) = terms;
    (*qtimgP) = qtimg;

    return;
    
    /* INTERNAL IMPLEMENTATIONS */

    char *describe_term_3x3(int32_t q)
      { 
         assert(NW == 3);

        /* The first 5 terms are the squares of all window ops except the average, condensed. */
        /* The next 28 terms are the products of distinct window ops except the average. */
        /* The last term is the constant 1.0. */

        if ((! cross_terms) && (q > 5)) { return NULL; }
        
        switch(q)
          { 
            case 0:
              return txtcat("1", "");
            case 1:
              return txtcat("P1^2+P2^2", ""); /* Gradient modulus squared. */
            case 2:
              return txtcat("P3^2+P4^2", ""); /* Second derivs squared. */
            case 3:
              return txtcat("P5^2", ""); /* Mixed derivative squared. */
            case 4:
              return txtcat("P6^2", ""); /* Checker squared. */
            case 5:
              return txtcat("P7^2+P8^2", ""); /* Saddle and co-saddle squared. */
            case 6:
              return txtcat("P2*P1", "");
            case 7:
              return txtcat("P3*P1", "");
            case 8:
              return txtcat("P3*P2", "");
            case 9:
              return txtcat("P4*P1", "");
            case 10:
              return txtcat("P4*P2", "");
            case 11:
              return txtcat("P4*P3", "");
            case 12:
              return txtcat("P5*P1", "");
            case 13:
              return txtcat("P5*P2", "");
            case 14:
              return txtcat("P5*P3", "");
            case 15:
              return txtcat("P5*P4", "");
            case 16:
              return txtcat("P6*P1", "");
            case 17:
              return txtcat("P6*P2", "");
            case 18:
              return txtcat("P6*P3", "");
            case 19:
              return txtcat("P6*P4", "");
            case 20:
              return txtcat("P6*P5", "");
            case 21:
              return txtcat("P7*P1", "");
            case 22:
              return txtcat("P7*P2", "");
            case 23:
              return txtcat("P7*P3", "");
            case 24:
              return txtcat("P7*P4", "");
            case 25:
              return txtcat("P7*P5", "");
            case 26:
              return txtcat("P7*P6", "");
            case 27:
              return txtcat("P8*P1", "");
            case 28:
              return txtcat("P8*P2", "");
            case 29:
              return txtcat("P8*P3", "");
            case 30:
              return txtcat("P8*P4", "");
            case 31:
              return txtcat("P8*P5", "");
            case 32:
              return txtcat("P8*P6", "");
            case 33:
              return txtcat("P8*P7", "");

            default:
              return NULL;
          }
      }
        
    double compute_term_3x3(int32_t ns, double vr[], int32_t q)
      {
         assert(ns == NS);
         assert(NW == 3);
         
         if ((! cross_terms) && (q > 5)) { assert(FALSE); }
        
         switch(q)
          { 
            case 0:
              return 1.0;
            case 1:
              return vr[1]*vr[1] + vr[2]*vr[2]; /* Gradient modulus squared. */
            case 2: 
              return vr[3]*vr[3] + vr[4]*vr[4]; /* Second derivs squared. */
            case 3:
              return vr[5]*vr[5]; /* Mixed derivative squared. */
            case 4:
              return vr[6]*vr[6]; /* Checker squared. */
            case 5:
              return vr[7]*vr[7] + vr[8]*vr[8]; /* Saddle and co-saddle squared. */
            case 6:
              return vr[2]*vr[1];
            case 7:
              return vr[3]*vr[1];
            case 8:
              return vr[3]*vr[2];
            case 9:
              return vr[4]*vr[1];
            case 10:
              return vr[4]*vr[2];
            case 11:
              return vr[4]*vr[3];
            case 12:
              return vr[5]*vr[1];
            case 13:
              return vr[5]*vr[2];
            case 14:
              return vr[5]*vr[3];
            case 15:
              return vr[5]*vr[4];
            case 16:
              return vr[6]*vr[1];
            case 17:
              return vr[6]*vr[2];
            case 18:
              return vr[6]*vr[3];
            case 19:
              return vr[6]*vr[4];
            case 20:
              return vr[6]*vr[5];
            case 21:
              return vr[7]*vr[1];
            case 22:
              return vr[7]*vr[2];
            case 23:
              return vr[7]*vr[3];
            case 24:
              return vr[7]*vr[4];
            case 25:
              return vr[7]*vr[5];
            case 26:
              return vr[7]*vr[6];
            case 27:
              return vr[8]*vr[1];
            case 28:
              return vr[8]*vr[2];
            case 29:
              return vr[8]*vr[3];
            case 30:
              return vr[8]*vr[4];
            case 31:
              return vr[8]*vr[5];
            case 32:
              return vr[8]*vr[6];
            case 33:
              return vr[8]*vr[7];

            default:
              assert(FALSE);
          }
      }
   }
