/* Last edited on 2007-10-28 19:44:42 by stolfi */

#include <cmp.h>

/* FINDING SEVERAL ALIGNMENTS AT THE SAME TIME */
  

vec_typedef(alignment_vec_t,alignment_vec,alignment_t);
  /* A list of candidate {alignment_t}s. */

alignment_vec_t find_alignment_mscale
  ( image_vec_t *iv, 
    int mult, 
    r2_t *dmax, 
    r2_t *rd
  );
  /* Computes up to {mult} distinct good alignments of the images in
    the list {iv} that locally minimize their mismatch.
    
    In each alignment, the displacement {d[i]} to be applied to the
    original reference point {iv.el[i].P} will be at most {±dmax.c[j]}
    along each axis {j}. The parameter {rd} defines the size of the
    relevant neighborhood in each image. */

alignment_vec_t enlarge_refine_and_prune_cands
  ( alignment_vec_t *av,
    image_vec_t *iv_r, 
    image_vec_t *iv_o, 
    int mult, 
    double step,
    r2_t *rd
  );
  /* Takes a list {av} of candidate alignments for a 
    set {img_r} of reduced-scale images, and converts them
    to candidate alignments for the original-scale
    images {img_o}. Then adjusts each alignment to provide
    a better match between the images of {iv_o}.
    Then sorts the resulting alignments by decreasing 
    quality, and retains only the best {mult} candidates.
    
    The given candidate alignments are {av.el[0..nc-1]} 
    where {nc = av.nel}.  The images are {iv_r.el[0..ni-1]}
    (reduced) and {iv_o.el[0..ni-1]} (original), where
    {ni = iv_r.nel = iv_o.nel}. 
    
    The reference point of each image in each alignment is adjusted
    by at most {step} in each direction. The image mismatch
    is computed with a windoing function of nominal radius {rd}.
    In any case, the reference point will not be allowed to 
    move outside the image's bounding box. */
 
void write_alignment_list(FILE *wr, alignment_vec_t *av, image_vec_t *iv);
  /* Writes to {wr} the list of alignments {av}, which are assumed to 
    refer to the list of images {iv}.  Each alignment is printed as
    its index in {av}, followed by two real coordinates for each
    image, all separated by spaces. */

int num_rel_shifts(int ni, int dmax, int sum);
  /* Computes the number of integer vectors in {[-dmax..+dmax]^ni}
    which add to {sum}. */
    
void generate_all_alignments(int ni, int dmax, double step, int nal, alignment_t *al);
  /* Generate all relative alignments for {ni} images whose entries 
    range in {[-dmax .. +dmax]} times {step}.
    Expects {nal} to be the exact number of such alignments. 
    Returns the alignments in {al[0..*nal-1]}. */

void recenter_displacements(int ni, double *dcol, double *drow);
  /* Shifts all displacements so that their sum is zero.  */

void adjust_alignment
  ( int ni, 
    image_vec_t *iml, 
    double dmax, 
    double epsilon, 
    double radius, 
    alignment_t *a
  );
  /* Adjusts the shifts {a->el[i]} for each image {i}, so
    as to minimize the mismatch of their neighborhoods. The maximum
    adjustment is {dmax} and the nominal precision required is {epsilon}, in
    each direction. The neighborhoods are defined by a Gaussian weight function of
    characteristic radius {radius}. A displacement of zero means to
    align the center of the image at the origin. The displacements are
    adjusted so as to add to zero.  */

double compare_cands(int ni, alignment_t *al, alignment_t *bl);
  /* Returns the maximum difference between the shifts of {al} and {bl}.  */
  
int compare_mismatches(interval_t *a, interval_t *b); 
  /* Returns +1 if {a} looks smaller than {b}, -1 if {a} looks larger than {b}, 
    and 0 if they look equal or incomparable. */

/* IMPLEMENTATIONS */

alignment_vec_t find_alignments_mscale
  ( image_vec_t *iv, 
    int mult, 
    r2_t *dmax, 
    r2_t *rd
  )
  {
    fprintf(stderr, "enter find_alignments_mscale\n");
    int ni = iv->nel;
    
    alignment_vec_t av;
    int ncan;
    int j, k;
    if ((dmax.c[0] < 0.5) && (dmax.c[1] < 0.5))
      { int nall = (int)(floor(maxsols+0.5));
        fprintf(stderr, "  computing all %d solutions...\n", nall);
        av = generate_all_alignments(ni, idmax, (r2_t){{ step, step }});
        ncan = nall;
      }
    else 
      { fprintf(stderr, "  recursion to smaller scale...\n");
        /* Generate {4*mult} solutions for coarser scale, then refine and select. */
        image_vec_t *iv_r = reduce_images(iv);
        /* Compute 4*mult solutions for coarse images */
        r2_t dmax_r;  r2_scale(0.5, &dmax, &dmax_r);
        r2_t rd_r; r2_scale(0.5, &rd, &rd_r);
        av = find_alignments_mscale(&iv_r, 4*mult, &dmax_r, step, &rd_r);
        enlarge_refine_and_prune_cands(&av, &iv_r, iv, mult, step, rd);
        free_images(&iv_r);
      }

    fprintf(stderr, "exit find_alignments_mscale\n");
    return av;
  }


alignment_vec_t enlarge_refine_and_prune_cands
  ( alignment_vec_t *av,
    image_vec_t *iv_r, 
    image_vec_t *iv_o, 
    int mult, 
    double step,
    r2_t *rd
  )
  {
    fprintf(stderr, "  rescaling and refinement...\n");
    int ni = iv->nel;
    for (k = 0; k < av->nel; k++) 
      { alignment_t *alg = &(av.el[k]);
        /* Map alignment {alg} from images {iv_r} to images {iv_o}: */
        int i;
        for (i=0; i<ni; i++)
          { image_t *img_r = &(iv_r->el[i]);
            image_t *img_o = &(iv_o->el[i]);
            map_point(&(alg->p[i]), img_r, img_o);
          }
        /* Refine {alg} by moving each point at most {step} from current pos: */
        adjust_alignment(alg, dmax, step, rd, iv);
      }

    /* Select {mult} best canutions: */
    fprintf(stderr, "  sorting and pruning...\n");
    sort_and_prune_alignments(av, step, mult);
  }


int num_rel_shifts(int ni, int dmax, int sum)
  {
    fprintf(stderr, "(%d)", ni);
    if (abs(sum) > ni*dmax)
      { return 0; }
    else if (ni <= 1) 
      { return 1; }
    else if (ni == 2) 
      { return 2*dmax + 1; }
    else if (ni == 3) 
      { return 3*dmax*(dmax + 1) + 1; }
    else if (ni == 4) 
      { return 8*dmax*(dmax*(2*dmax + 3) + 1)/3 + 1; }
    else if (ni == 5) 
      { /* There must be a better way... */
        int i, s = 0;
        for (i=-dmax; i <= dmax; i++) { s += num_rel_shifts(ni-1,dmax,sum-i); }
        return s;
      }
    else
      { demand(FALSE, "too many images");
        return 0;
      }
  }

void sort_and_prune_alignments(int ni, int *ncanP, alignment_t *can, double epsilon, int max)
  {
    /* !!! must eliminate duplicates !!! */
    int ncan = *ncanP;
    int nok, i;
    /* Select the best {max} alignments: */
    nok = 1;
    for (i = 1; i < ncan; i++)
      { int j;
        alignment_t t = can[i];
        j = nok-1;
        while ((j >= 0) && (compare_cands(ni, &(can[j]), &t) > epsilon/2)) { j--; }
        if (j < 0)
          { j = nok;
            while ((j > 0) && (compare_mismatches(&(can[j-1].tvar), &(t.tvar)) > 0))
              { if (j < max) { can[j] = can[j-1]; }
                j--;
              }
            can[j] = t; if (nok < max) { nok++; }
          }
      }
    (*ncanP) = nok;
  }

double compare_cands(int ni, alignment_t *al, alignment_t *bl)
  {
    int i;
    double dmax = 0;
    for (i = 0; i < ni; i++)
      { double dc = fabs(al->d[i].c[0] - bl->d[i].c[0]);
        double dr = fabs(al->d[i].c[1] - bl->d[i].c[1]);
        if (dc > dmax) { dmax = dc; }
        if (dr > dmax) { dmax = dr; }
      }
    return dmax; 
  }

int compare_mismatches(interval_t *a, interval_t *b)
  {
    double amin = 0.75*a->end[0] + 0.25*a->end[1];
    double bmin = 0.75*b->end[0] + 0.25*b->end[1];
    return cmp_double(&amin, &bmin); 
  }

alignment_t *alloc_alignment_set(int ncan)
  {
    int i;
    alignment_t *als = (alignment_t *)malloc(ncan*sizeof(alignment_t));
    if (als == NULL) { demand(FALSE, "out of memory for alignments"); }
    for (i = 0; i < ncan; i++) 
      { alignment_t *alsi = &(als[i]);
        alsi->d = NULL; alsi->mism = NULL; 
        alsi->tvar = (interval_t){0.00, 0.25};
      }
    return als;
  }

void generate_all_alignments(int ni, int dmax, double step, int nal, alignment_t *al)
  {
    int sc, sr;
    int *dc = (int *)malloc(ni*sizeof(int));
    int *dr = (int *)malloc(ni*sizeof(int));
    int k;
    int carry; /* boolean */
    int nfound = 0;

    k = ni; sc = 0; sr = 0; carry = 0;
    while(1)
      { 
        /* Here {sc} is the sum of {dc[k..ni-1]}, ditto for {sr,dr}. */
        /* Generate the smallest solution with current {dc,dr[k..nfound-1]}: */
        int j;
        while (k > 0)
          { /* Move down and set digit to minimum value: */
            k--;
            dc[k] = imax(-dmax, -(k*dmax + sc));
            dr[k] = imax(-dmax, -(k*dmax + sr));
            sc += dc[k]; sr += dr[k];
          }     

        /* Store this solution: */
        if (nfound >= nal) 
          { demand(FALSE, "num of alignments underestimated"); }
        if ((sc != 0) || (sr != 0)) 
          { demand(FALSE, "displacements don't add to zero"); }
        { alignment_t a = alloc_alignment(ni);
          fprintf(stderr, "\n");
          for (j = 0; j < ni; j++)
            { a.d[j] = (r2_t){{ step*dc[j], step*dr[j] }};
              a.mism[j] = (interval_t){0.00, 0.25};
              debug_displacement("  ", j, a.d[j], a.mism[j], "\n");
            }
          a.tvar = (interval_t){0.00, 0.25};
          al[nfound] = a; nfound++;
        }

        /* Find first position {k} that can be incremented, and increment it: */
        carry = 1;
        while (carry)
          { /* Now {sc,sr} is sum of {dc,dr[k..ni-1]} */
            if (k >= ni)
              { /* Exhausted all possibilities, stop: */
                if (nfound != nal) { demand(FALSE, "number of alignments doesn't check"); }
                free(dc); free(dr);
                return;
              }
            sc -= dc[k]; sr -= dr[k]; 
            if (dr[k] < imin(dmax, k*dmax - sr))
              { /* Bump {dr[k]}, reset carry: */
                dr[k]++;
                carry = 0; 
              }
            else if (dc[k] < imin(dmax, k*dmax - sc))
              { /* Bump {dc[k]}, reset {dr[k]}, reset carry: */
                dc[k]++; dr[k] = imax(-dmax, -(k*dmax+sr)); 
                carry = 0; 
              }
            else 
              { /* This position can't be incremented, carry to the next one: */
                k++; 
              }
            /* Now {sc,sr} is sum of {dc,dr[k+1-carry..ni-1]} */
          }
        sc += dc[k]; sr += dr[k];
      }
  }

void adjust_alignment
  ( int ni, 
    image_vec_t *iml, 
    double dmax, 
    double epsilon, 
    double radius, 
    alignment_t *a
  )
  {
    if (ni >= 2)
      { 
        int arad = ceil(3*radius); /* Weights are negligible beyond this. */
        int nwt = 2*arad+1;
        float_image_t *avg = allocate_image(nwt, nwt, max_chns(iml));
        double *gwt = gaussian_distr(nwt, radius);
        double s1 = ((double)(ni-1))/((double)ni);
        double corr = s1*s1; /* Variance correction factor. */
        image_vec_t *impi;
        double step;
        int i;

        /* Loop until convergence: */
        i = 0; impi = iml;
        step = 0.5*(dmax + epsilon);
        while(1)
          { int krowi, kcoli;
            interval_t prev_best; 
            /* Adjust displacements so that they add to zero: */
            recenter_displacements(ni, al->dcol, al->drow);
            /* Compute the average of all images excluding image {i}: */
            fprintf(stderr, "%s: averaging images except %d...", PROG_NAME, i);
            average_images(ni, iml, i, avg, al->dcol, al->drow);
            if (avg->cols < 20) { debug_image("avg", "", avg); }
            /* Compare image {i} against average, in all possible displacements: */
            fprintf(stderr, "%s: computing current mismatch of image %d...", PROG_NAME, i);
            al->mism[i] = compute_mismatch(impi->im, al->d[i].c[0], al->d[i].c[1], avg, corr, nwt, gwt);
            debug_displacement("INIT", i, al->d[i].c[0], al->d[i].c[1], al->mism[i], "\n");
            prev_best = al->mism[i];
            for (krowi = -1; krowi <= +1; krowi++)
              { double drowi = al->d[i].c[1] + step*krowi;
                for (kcoli = -1; kcoli <= +1; kcoli++)
                  { double dcoli = al->d[i].c[0] + step*kcoli;
                    interval_t s = compute_mismatch(impi->im, dcoli, drowi, avg, corr, nwt, gwt);
                    if (compare_mismatches(&s, &(al->mism[i])) < 0) 
                      { debug_displacement("good", i, dcoli, drowi, s, "\n");
                        al->mism[i] = s; al->d[i].c[0] = dcoli; al->d[i].c[1] = drowi;
                      }
                  }
              }
            if (compare_mismatches(&(al->mism[i]), &prev_best) < 0) 
              { debug_displacement("BEST", i, al->d[i].c[0], al->d[i].c[1], al->mism[i], "\n");
              }
            impi = impi->next; i++;
            if (i >= ni) 
              { i = 0; impi = iml; 
                if (step <= 1.00001*epsilon) 
                  { recompute_pixel_variance(ni, iml, al, avg, nwt, gwt); return; } 
                step /= 2;
              }
          }
      }
  }

void recompute_pixel_variance(int ni, image_vec_t *iml, alignment_t *al, float_image_t *avg, int nwt, double *wt)
  {
    int i = 0;
    image_vec_t *impi = iml;
    double qt = 1.0/((double)ni);
    al->tvar = (interval_t){0, 0};
    average_images(ni, iml, ni, avg, al->dcol, al->drow);
    for (i = 0, impi=iml; i < ni; i++, impi = impi->next)
      { interval_t s = compute_mismatch(impi->im, al->d[i].c[0], al->d[i].c[1], avg, 1.0, nwt, wt);
        al->mism[i] = s;
        al->tvar.end[0] += qt*s.end[0];
        al->tvar.end[1] += qt*s.end[1];
      }
  }

void recenter_displacements(int ni, double *dcol, double *drow)
  {
    double tcol = 0, trow = 0;
    int i;
    for (i = 0; i < ni; i++) { tcol += d[i].c[0]; trow += d[i].c[1]; }
    tcol /= (double)ni;
    trow /= (double)ni;
    for (i = 0; i < ni; i++) { d[i].c[0] -= tcol; d[i].c[1] -= trow; }
  }
