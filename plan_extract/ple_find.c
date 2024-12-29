/* See ple_find.h */
/* Last edited on 2024-12-21 14:18:51 by stolfi */


#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <i2.h>
#include <irange.h>
#include <affirm.h>

#include <ple_find.h>

/* INTERNAL PROTOTYPES */

void ple_eval_plans_for_ssize_esize
  ( irange_t TBox[], 
    bool_t paddable[], 
    i2_t *SSize, 
    i2_t *ESize, 
    i2_t *EStep, 
    double maxScale, 
    double maxStretch, 
    double maxPad, 
    double maxCut, 
    plan_t *pl_best, 
    int debug
  );
  /* Similar to {ple_eval_plans_for_esize}, but considers only plans where the 
    image {S} has dimensions {SSize}. */
    
void ple_eval_plans_for_psize_ssize_esize
  ( irange_t TBox[], 
    bool_t paddable[], 
    i2_t *PSize, 
    i2_t *SSize, 
    i2_t *ESize, 
    i2_t *EStep, 
    double maxScale, 
    double maxStretch, 
    double maxPad, 
    double maxCut, 
    plan_t *pl_best, 
    int debug
  );
  /* Similar to {ple_eval_plans_for_ssize_esize}, but considers only plans where the 
    image {P} has dimensions {PSize}. */
    
void ple_compute_plan_score
  ( irange_t TBox[], 
    bool_t paddable[], 
    plan_t *pl, 
    double maxScale, 
    double maxStretch, 
    double maxPad, 
    double maxCut, 
    int debug
  );
  /* (Re)computes the score of plan {pl}, saves it in {pl->score}.
    The domain {pl->CBox} of the {C} image must be contained in {TBox}.
    Assumes that {T} can be padded along edge {e} iff {paddable[e]} is true. */

double ple_cut_factor(irange_t T[], irange_t C[]);
  /* The amount of box {T} that is not contained in box {C},
    as a fraction of {T}'s extent, max-ed over both axes. 
    Assumes that {C} is contained in {T}. */

double ple_pad_factor(irange_t C[], irange_t P[]);
  /* The amount of box {P} that is not contained in box {C},
    as a fraction of {P}'s extent, max-ed over both axes. 
    Assumes that {P} is relative to {C} and contains it. */

void ple_eval_plans_for_esize
  ( irange_t TBox[], 
    bool_t paddable[], 
    i2_t *ESize, 
    i2_t *EStep,
    double maxScale, 
    double maxStretch, 
    double maxPad, 
    double maxCut, 
    plan_t *pl_best, 
    int debug
  )
  {
    demand(maxPad < 1.0, "invalid {maxPad}");
    
    bool_t verbose = (debug > 0);
    if (verbose) 
      { fprintf(stderr, "+ ple_eval_plans_for_esize\n"); 
        ple_debug_size(stderr, "  ESize = ", ESize, "\n");
        ple_debug_size(stderr, "  EStep = ", EStep, "\n");
        fprintf(stderr, "\n");
      }
    
    /* Exclude some degenerate cases: */
    if ((ESize->c[0] <= 0) || (ESize->c[1] <= 0)) { return; }

    /* Size of {T} sub-image: */
    i2_t TSize = ple_box_size(TBox);
    
    /* Find the min and max dimensions {minPSize,maxPSize} of {P}, for any {S}. */
    i2_t minPSize;       /* Minimum dimensions of {P}. */
    i2_t maxPSize;       /* Maximum dimensions of {P}. */
    int ax;
    for (ax = 0; ax < 2; ax++)
      { /* The minimum size of {P} is determined by {maxScale} and {maxCut}: */
        int min1 = (int)floor(ESize->c[ax]/maxScale); 
        int min2 = (int)floor((1-maxCut)*TSize.c[ax]);
        minPSize.c[ax] = (min1 > min2 ? min1 : min2);
        /* Get the two edges perpendicular to {ax}: */
        edge_t eLo = (ax == 0 ? edge_L : edge_B);
        edge_t eHi = (ax == 0 ? edge_R : edge_T);
        /* Determine how much {T} can be padded along axis {ax}: */
        double maxPRelSize = (paddable[eLo] || paddable[eHi] ? 1/(1 - maxPad) : 1); 
        /* Max extent of {P} is that of {T} times {maxPRelSize}: */
        maxPSize.c[ax] = (int)ceil(maxPRelSize*TSize.c[ax]);
        /* Check against excessive size: */
        if (maxPSize.c[ax] > ple_MAX_SIZE) { maxPSize.c[ax] = ple_MAX_SIZE; }
      }
    if (verbose) 
      { ple_debug_size(stderr, "  maxPSize = ", &maxPSize, "\n");
        ple_debug_size(stderr, "  minPSize = ", &minPSize, "\n");
      }

    /* Determine the min scale factor {minScale} for the {P->S} scaling: */
    double minScale = 0.0;      /* Minimum linear factor for {P->S} scaling. */
    for (ax = 0; ax < 2; ax++)
      { /* Get scale factor {r} assuming one tile and the largest possible {P}: */
        double r = 0.999999*((double)ESize->c[ax])/((double)maxPSize.c[ax]);
        if (r > minScale) { minScale = r; }
      }
    if (verbose) 
      { fprintf(stderr, "  minScale = %7.5f\n", minScale); }

    /* Get min and max tile count {ntiles[0..1]} and tile step {EStep}: */
    irange_t ntiles[2];
    i2_t minSSize, maxSSize; /* Min and max size of {S} image. */
    for (ax = 0; ax < 2; ax++)
      { /* The minimum size is defined by {minPSize} and {minScale}: */
        minSSize.c[ax] = (int)floor(minScale*minPSize.c[ax]);
        /* The maximum size is defined by {maxPSize} and {maxScale}: */
        maxSSize.c[ax] = (int)ceil(maxScale*maxPSize.c[ax]);
        /* Provide sensible limits for degenerate cases: */
        if (minSSize.c[ax] < ESize->c[ax]) { minSSize.c[ax] = ESize->c[ax]; }
        if (maxSSize.c[ax] > ple_MAX_SIZE) { maxSSize.c[ax] = ple_MAX_SIZE; }
        /* The {S} size increment is twice {step}: */
        int incr = 2*EStep->c[ax];
        /* Get the min tile count {ntiles[ax].end[0]} from {minSSize}, rounding up: */
        int dmin = minSSize.c[ax] - ESize->c[ax];
        assert(dmin >= 0);
        ntiles[ax].end[0] = 1 + 2*((dmin + incr - 1)/incr);
        /* Get the max tile count {ntiles[ax].end[1]} from {maxSSize}, rounding down: */
        int dmax = maxSSize.c[ax] - ESize->c[ax];
        ntiles[ax].end[1] = (dmax < 0 ? -1 : 1 + 2*(dmax/incr));
      }
    if (verbose) 
      { ple_debug_size(stderr, "  minSSize = ", &minSSize, "\n");
        ple_debug_size(stderr, "  maxSSize = ", &maxSSize, "\n");
        ple_debug_ixrange(stderr, "  nX = ", &(ntiles[0]), "\n");
        ple_debug_ixrange(stderr, "  nY = ", &(ntiles[1]), "\n");
      }

    /* Consider all odd tile counts {ntx,nty} allowed by {ntiles[0..1]}: */
    int ntx, nty;
    for (ntx = ntiles[0].end[0]; ntx <= ntiles[0].end[1]; ntx += 2)
      { for (nty = ntiles[1].end[0]; nty <= ntiles[1].end[1]; nty += 2)
          { /* Compute size of {S} from {ntx, nty}: */
            i2_t SSize;
            SSize.c[0] = ESize->c[0] + (ntx - 1)*EStep->c[0]; 
            SSize.c[1] = ESize->c[1] + (nty - 1)*EStep->c[1]; 
            /* Evaluate the plans with this {SSize,ESize}: */
            ple_eval_plans_for_ssize_esize
              ( TBox, paddable, &SSize, ESize, EStep,
                maxScale, maxStretch, maxPad, maxCut, 
                pl_best, debug - 10
              );
          }
      }

    if (verbose) 
      { fprintf(stderr, "\n");
        fprintf(stderr, "- ple_eval_plans_for_esize\n");
      }
  }

void ple_eval_plans_for_ssize_esize
  ( irange_t TBox[], 
    bool_t paddable[], 
    i2_t *SSize, 
    i2_t *ESize, 
    i2_t *EStep, 
    double maxScale, 
    double maxStretch, 
    double maxPad, 
    double maxCut, 
    plan_t *pl_best, 
    int debug
  )
  {
    bool_t verbose = (debug > 0);
    if (verbose) 
      { fprintf(stderr, "  + ple_eval_plans_for_ssize_esize\n"); 
        ple_debug_size(stderr, "    SSize = ", SSize, "\n");
        ple_debug_size(stderr, "    ESize = ", ESize, "\n");
        ple_debug_size(stderr, "    EStep = ", EStep, "\n");
        fprintf(stderr, "\n");
      }
    
    /* Exclude some degenerate cases: */
    if ((ESize->c[0] <= 0) || (ESize->c[1] <= 0)) { return; }
    if ((SSize->c[0] <= 0) || (SSize->c[1] <= 0)) { return; }

    /* Size of {T} sub-image: */
    i2_t TSize = ple_box_size(TBox);
    
    /* Find the longest and shortest axes of {SSize}: */
    int axMax = (SSize->c[0] >= SSize->c[1] ? 0 : 1);
    int axMin = 1 - axMax;

    /* Find the min and max dimensions {minPSize,maxPSize} of {P} given {SSize}. */
    i2_t minPSize;   /* Minimum dimensions of {P}. */
    i2_t maxPSize;   /* Maximum dimensions of {P}. */
    int ax;
    for (ax = 0; ax < 2; ax++)
      { /* The minimum size of {P} is determined by {maxScale} and {maxCut}: */
        int min1 = (int)floor(SSize->c[ax]/maxScale); 
        int min2 = (int)floor((1-maxCut)*TSize.c[ax]);
        minPSize.c[ax] = (min1 > min2 ? min1 : min2);
        /* Get the two edges perpendicular to {ax}: */
        edge_t eLo = (ax == 0 ? edge_L : edge_B);
        edge_t eHi = (ax == 0 ? edge_R : edge_T);
        /* Determine how much {T} can be padded along axis {ax}: */
        double maxPRelSize = (paddable[eLo] || paddable[eHi] ? 1/(1 - maxPad) : 1); 
        /* Max extent of {P} cannot exceed that of {T} times {maxPRelSize}: */
        maxPSize.c[ax] = (int)ceil(maxPRelSize*TSize.c[ax]);
        /* Provide sensible limits for degenerate cases: */
        if (minPSize.c[ax] < 1) { minPSize.c[ax] = 1; }
        if (maxPSize.c[ax] > ple_MAX_SIZE) { maxPSize.c[ax] = ple_MAX_SIZE; }
      }

    /* Consider all alternatives for the scaled image size {PSize}: */
    int tad = 1; /* Fudge longest side by this much around ideal value. */
    i2_t PSize = minPSize;
    double aspect = ((double)SSize->c[axMax])/((double)SSize->c[axMin]);
    while ((PSize.c[0] <= maxPSize.c[0]) && (PSize.c[1] <= maxPSize.c[1]))
      { /* Evaluate the plan: */
        ple_eval_plans_for_psize_ssize_esize
          ( TBox, paddable, &PSize, SSize, ESize, EStep,
            maxScale, maxStretch, maxPad, maxCut, 
            pl_best, debug - 10
          );
        /* Compute upper limits for {PSize.c[axMax]} given {PSize.c[axMin]}: */
        int max1 = maxPSize.c[axMax];
        int max2 = (int)ceil(aspect*PSize.c[axMin]) + tad;
        /* Increment {PSize} along longest axis, if possible: */
        PSize.c[axMax]++;
        if ((PSize.c[axMax] > max1) || (PSize.c[axMax] > max2))
          { /* Cannot increment {PSize.c[axMax]} any more. */
            /* Increment {PSize.c[axMin]} instead: */
            PSize.c[axMin]++;
            /* Compute lower limits for {PSize.c[axMax]} given {PSize.c[axMin]}: */
            int min1 = minPSize.c[axMax];
            int min2 = (int)floor(aspect*PSize.c[axMin]) - tad;
            PSize.c[axMax] = (min1 > min2 ? min1 : min2);
          }
      }

    if (verbose) 
      { fprintf(stderr, "\n");
        fprintf(stderr, "  - ple_eval_plans_for_ssize_esize\n");
      }
  }

void ple_eval_plans_for_psize_ssize_esize
  ( irange_t TBox[], 
    bool_t paddable[], 
    i2_t *PSize, 
    i2_t *SSize, 
    i2_t *ESize, 
    i2_t *EStep, 
    double maxScale, 
    double maxStretch, 
    double maxPad, 
    double maxCut, 
    plan_t *pl_best, 
    int debug
  )
  {
    plan_t pl;
    pl.SSize = *SSize;
    pl.ESize = *ESize;
    pl.EStep = *EStep;
    
    bool_t verbose = (debug > 0);
    if (verbose) 
      { fprintf(stderr, "  + ple_eval_plans_for_psize_ssize_esize\n"); 
        ple_debug_size(stderr, "    PSize = ", PSize, "\n");
        ple_debug_size(stderr, "    SSize = ", SSize, "\n");
        ple_debug_size(stderr, "    ESize = ", ESize, "\n");
        ple_debug_size(stderr, "    EStep = ", EStep, "\n");
        fprintf(stderr, "\n");
      }

    /* Size of {T} sub-image: */
    i2_t TSize = ple_box_size(TBox);
    
    /* Compute {pl.PBox}, {pl.CBox}, and {pl.EPos}: */
    int ax;
    for (ax = 0; ax < 2; ax++) 
      { /* Define {pl.CBox[ax]}, {pl.PBox[ax]}: */
        if (PSize->c[ax] <= TSize.c[ax])
          { /* The best {pl.CBox} is centered in {TBox} along {ax}: */
            pl.CBox[ax].end[0] = TBox[ax].end[0] + (TSize.c[ax] - PSize->c[ax])/2;
            pl.CBox[ax].end[1] = pl.CBox[ax].end[0] + PSize->c[ax];
            /* No padding needed: */
            pl.PBox[ax].end[0] = 0;
            pl.PBox[ax].end[1] = pl.PBox[ax].end[0] + PSize->c[ax];
          }
        else 
          { /* The best {pl.CBox} is the whole of {TBox} along {ax}: */
            pl.CBox[ax] = TBox[ax];
            /* Get the two edges perpendicular to {ax}: */
            edge_t eLo = (ax == 0 ? edge_L : edge_B);
            edge_t eHi = (ax == 0 ? edge_R : edge_T);
            /* Define the low edge of {pl.PBox} relative to {pl.CBox}: */
            if ((paddable[eLo]) && (paddable[eHi]))
              { /* Pad equally on both sides (best choice): */
                pl.PBox[ax].end[0] = - (PSize->c[ax] - TSize.c[ax])/2;
              }
            else if (paddable[eLo])
              { /* Pad the low side only: */
                pl.PBox[ax].end[0] = - (PSize->c[ax] - TSize.c[ax]);
              }
            else if (paddable[eHi])
              { /* Pad the high side only: */
                pl.PBox[ax].end[0] = 0;
              }
            else
              { /* The domain {TBox} cannot be padded to {PSize}: */
                return; 
              }
            /* Define the high edge of {pl.PBox}: */
            pl.PBox[ax].end[1] = pl.PBox[ax].end[0] + PSize->c[ax];
          }
      }
                
    /* Consider each axis in turn: */
    for (ax = 0; ax < 2; ax++) 
      { /* Make sure that {S} is large enough to extract at least one tile: */
        if (ESize->c[ax] > SSize->c[ax]) { return; }
        /* Define the position {pl.EPos[ax]} of {E[0,0]} rel to {S}: */
        pl.EPos.c[ax] = (SSize->c[ax] - ESize->c[ax])/2;
        /* Make sure that {S} precisely contains an odd number {2*N+1} of tiles: */
        int W = SSize->c[ax] - ESize->c[ax];
        int N = W / (2*EStep->c[ax]);
        if (N * (2*EStep->c[ax]) != W) { return; }
        /* Compute the range for {v} that covers {SSize} along {ax}: */
        irange_t v = (irange_t){{ -N, +N }}; /* Range that covers the full extent: */
        assert ((v.end[0] <= 0) && (v.end[1] >= 0));
        pl.tr[ax] = v;
      }

    /* Evaluate and update: */
    ple_compute_plan_score
      ( TBox, paddable, &pl, maxScale, maxStretch, maxPad, maxCut, debug - 10 );
    if (pl.score > pl_best->score) { (*pl_best) = pl; }

    if (verbose) 
      { fprintf(stderr, "\n");
        fprintf(stderr, "  - ple_eval_plans_for_psize_ssize_esize\n");
      }
  }

void ple_compute_plan_score
  ( irange_t TBox[], 
    bool_t paddable[], 
    plan_t *pl, 
    double maxScale, 
    double maxStretch, 
    double maxPad, 
    double maxCut, 
    int debug
  )
  {
    /* Size of {T,C,P} sub-images: */
    i2_t TSize = ple_box_size(TBox);
    i2_t CSize = ple_box_size(pl->CBox);
    i2_t PSize = ple_box_size(pl->PBox);

    bool_t verbose = (debug > 0);
    if (verbose)
      { ple_debug_plan(stderr, "plan:\n", TBox, pl, "\n"); }

    /* Tile counts and tile-count factor: */
    int ntiles[2]; /* Number of tiles along each axis. */
    int ax;
    for (ax = 0; ax < 2; ax++) 
      { /* Compute number of tiles in each direction: */
        ntiles[ax] = (pl->tr[ax].end[1] - pl->tr[ax].end[0] + 1);
      }
    int nt = ntiles[0]*ntiles[1];
    double tileFactor = 1.0 - 0.5*(nt - 1); 
    if (verbose) { fprintf(stderr, "  tileFactor = %7.5f\n", tileFactor); }
    
    /* Compute partial scores: 1 for best possible, negative if constraint is violated. */

    /* Compute {scaleScore}: Does {S/P} exceed {maxScale}?  */
    double scaleScore = 1.0; /* Score for {maxScale} compliance. */
    for (ax = 0; ax < 2; ax++) 
      { int MDim = (int)ceil(maxScale*PSize.c[ax]); /* Max size for {S} along {ax}. */
        int SDim = pl->SSize.c[ax]; /* Actual size of {S} along {ax}. */
        double s = (SDim > MDim ? -1.0 : +1.0);
        if (s < scaleScore) { scaleScore = s; }
      }
    if (verbose) { fprintf(stderr, "  scaleScore = %7.5f\n", scaleScore); }

    /* Compute {shapeScore}: How close is the aspect of {S} to that of {P}? */
    double aspS = ((double)pl->SSize.c[0])/((double)pl->SSize.c[1]);
    double aspP = ((double)PSize.c[0])/((double)PSize.c[1]);
    double r = log(aspS/aspP);
    double shapeScore = 1 - (r*r)/(2*log(maxStretch));
    if (verbose) { fprintf(stderr, "  shapeScore = %7.5f\n", shapeScore); }

    /* Compute {padScore}: Is {maxPad} honored? How much padding is needed? */
    double padScore = 1.0;
    for (ax = 0; ax < 2; ax++) 
      { /* Compute the amount of padding required at each end of {P}: */
        int padLo = - pl->PBox[ax].end[0];
        affirm(padLo >= 0, "bug padLo");
        int padHi = pl->PBox[ax].end[1] - CSize.c[ax];
        affirm(padHi >= 0, "bug padHi");
        /* Compute actual padding ratio {padR} along this axis: */
        int padR = ((double)padLo + padHi)/((double)PSize.c[ax]);
        /* Get the edges {eLo,eHi} perpendicular to {ax}: */
        edge_t eLo = (ax == 0 ? edge_L : edge_B);
        edge_t eHi = (ax == 0 ? edge_R : edge_T);
        /* Compute pad score {s} for this axis: */
        double s;
        if (((! paddable[eLo]) && (padLo > 0)) || ((! paddable[eHi]) && (padHi > 0)))
          { /* Image does not allow padding as planned: */
            s = -1.0;
          }
        else if (padR <= 0.0)
          { /* Plan requires no padding: */
            s = +1.0;
          }
        else if (padR > maxPad)
          { /* Plan requires too much padding: */
            s = -1.0;
          }
        else
          { /* Compute the ratio {r} of planned padding to max padding: */
            double r = ((padR - 1)/maxPad);
            s = 1.0 - r*r;
            if (s < 0) { s = 0.0; }
          }
        /* Update {padScore}: */
        if (s < padScore) { padScore = s; }
      }
    if (verbose) { fprintf(stderr, "  padScore =   %7.5f\n", padScore); }

    /* Compute {cutScore}: How much of {TBox} is left uncovered? */
    double cutScore = 1.0;
    for (ax = 0; ax < 2; ax++) 
      { /* Compute the amount of cropping required at each end: */
        int cutLo = pl->CBox[ax].end[0] - TBox[ax].end[0];
        int cutHi = TBox[ax].end[1] - pl->CBox[ax].end[1];
        /* Compute the maximum amount of cropping allowed along {ax}: */
        int numCut = (int)ceil(maxCut * TSize.c[ax]); 
        /* Compute the cut score {s} for this axis: */
        double s;
        if (cutLo + cutHi > numCut)
          { /* Plan requires too much cropping: */
            s = -1.0;
          }
        else
          { /* Compute the ratio {r} of planned cropping to max cropping: */
            double r = ((double)cutLo + cutHi)/((double)numCut);
            s = 1.0 - r*r;
            if (s < 0) { s = 0.0; }
          }
        /* Update {cutScore}: */
        if (s < cutScore) { cutScore = s; }
      }
    if (verbose) { fprintf(stderr, "  cutScore =   %7.5f\n", cutScore); }

    /* Compute combined score = minimum of partial scores. */
    double score = 1.0;
    if (scaleScore < score) { score = scaleScore; }
    if (shapeScore < score) { score = shapeScore; }
    if (padScore < score) { score = padScore; }
    if (cutScore < score) { score = cutScore; }

    /* If positive, multiply by {tileFactor}: */
    if (score > 0) { score *= tileFactor; }
    if (verbose) { fprintf(stderr, "  score =      %7.5f\n", score); }
    
    /* Save score: */
    pl->score = score;
  }

void ple_debug_plan(FILE *wr, char *pre, irange_t TBox[], plan_t *pl, char *suf)
  {
    /* Debugging printout: */
    fprintf(wr, "%s:\n", pre); 
    ple_debug_box(wr,     "  TBox =   ", TBox, "\n");
    ple_debug_box(wr,     "  CBox =   ", pl->CBox, "");
    double cutFactor = ple_cut_factor(TBox, pl->CBox);
    fprintf(wr, "  cutFactor = %7.5f\n", cutFactor);
    ple_debug_box(wr,     "  PBox =   ", pl->PBox, "");
    double padFactor = ple_pad_factor(pl->CBox, pl->PBox);
    fprintf(wr, "  padFactor = %7.5f\n", padFactor);
    ple_debug_size(wr,    "  SSize =  ", &(pl->SSize), "\n");
    ple_debug_size(wr,    "  ESize =  ", &(pl->ESize), "\n");
    ple_debug_pos(wr,     "  EPos =   ", &(pl->EPos), "\n");
    ple_debug_pos(wr,     "  EStep =  ", &(pl->EStep), "\n");
    ple_debug_ixrange(wr, "  rrange = ", &(pl->tr[0]), "\n");
    ple_debug_ixrange(wr, "  srange = ", &(pl->tr[1]), "\n");
    fprintf(wr, "%s", suf); 
  }

void ple_debug_box(FILE *wr, char *pre, irange_t B[], char *suf) 
  { fprintf(wr, "%s", pre);
    fprintf(wr, "[%d..%d] × [%d..%d]", 
      B[0].end[0], B[0].end[1], B[1].end[0], B[1].end[1]
    );
    fprintf(wr, "%s", suf);
  }

void ple_debug_size(FILE *wr, char *pre, i2_t *sz, char *suf)
  { fprintf(wr, "%s", pre);
    fprintf(wr, "%d × %d", sz->c[0], sz->c[1]);
    fprintf(wr, "%s", suf);
  }

void ple_debug_pos(FILE *wr, char *pre, i2_t *pos, char *suf)
  { fprintf(wr, "%s", pre);
    fprintf(wr, "( %d %d )", pos->c[0], pos->c[1]);
    fprintf(wr, "%s", suf);
  }

void ple_debug_ixrange(FILE *wr, char *pre, irange_t *v, char *suf)
  { fprintf(wr, "%s", pre);
    fprintf(wr, "{ %d .. %d }", v->end[0], v->end[1]);
    fprintf(wr, "%s", suf);
  }

i2_t ple_box_size(irange_t B[])
  { return 
      (i2_t){{ 
          B[0].end[1] - B[0].end[0],
          B[1].end[1] - B[1].end[0]
        }};
  }

double ple_cut_factor(irange_t T[], irange_t C[])
  { int ax;
    double cutFactor = 0;
    for (ax = 0; ax < 2; ax++)
      { int TSize = T[ax].end[1] - T[ax].end[0];
        int CSize = C[ax].end[1] - C[ax].end[0];
        demand(C[ax].end[0] >= T[ax].end[0], "{C} starts before {T}");
        demand(C[ax].end[1] <= T[ax].end[1], "{C} stops after {T}");
        double f = ((double)(TSize - CSize))/((double)TSize);
        if (f > cutFactor) { cutFactor = f; }
      }
    return cutFactor;
  }

double ple_pad_factor(irange_t C[], irange_t P[])
  { int ax;
    double padFactor = 0;
    for (ax = 0; ax < 2; ax++)
      { int CSize = C[ax].end[1] - C[ax].end[0];
        int PSize = P[ax].end[1] - P[ax].end[0];
        demand(P[ax].end[0] <= 0, "{P} starts after {C}");
        demand(P[ax].end[1] >= CSize, "{P} stops before {C}");
        double f = ((double)(PSize - CSize))/((double)PSize);
        if (f > padFactor) { padFactor = f; }
      }
    return padFactor;
  }

bool_t ple_empty_box(irange_t B[])
  { return 
      ( B[0].end[0] >= B[0].end[1]) ||
      ( B[1].end[0] >= B[1].end[1]);
  }

