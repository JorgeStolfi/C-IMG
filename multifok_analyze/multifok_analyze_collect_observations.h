#ifndef multifok_analyze_collect_observations_H
#define multifok_analyze_collect_observations_H
/* Last edited on 2018-09-05 17:58:10 by stolfilocal */

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <float_image.h>

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
  );
  /* Collects observations to infer the optimal focus operator by quadratic regression.
  
    The input is a reference (totally focused) image {ref} and a set of
    {NF} frames {frame[0..NF-1]} with limited depth of focus, all
    monochromatic and with the same dimensions. The output is a set of
    {NP} observations for linear regression, each consisting of values of
    a dependent variable {e2[p]} and {NQ} independent variables
    {terms[p][0..NQ-1]}, for {p} in {0..NP-1}. The counts {NP} and {NQ}
    are determined by the procedure.
    
    More precisely, sweeps a sliding window {NW} by {NW} pixels over the
    frame images and the reference image, considering all placements
    such that the window fits entirely in the image. 
    
    For each such window position extracts the {NS=NW*NW}
    samples from that window of the reference image as a vector
    {rr[0..NS-1]}. If the deviation of those samples is at least
    {3*noise}, that window position is considered significant and will
    generate {NF} observations. Otherwise that window position is
    ignored.
    
    Let {ND} be the number of all significant window positions. The
    number of observations {NP} will be {ND*NF}.
    
    The following steps are performed for each significant window
    position. The vector {rr[0..NS-1]} is normalized to zero mean and
    unit norm.  Then, for each frame index {f} in {0..NF-1}, the procedure extracts
    the {NS} samples in that window of the {frame[f]} as a vector
    {fr[0..NS-1]}, and normalizes that vector too.
    
    The procedure then computes the quadratic discrepancy {E[p]},
    defined as the distance squared of the two
    vectors, and stores it in {(*e2P)[p]), where {p=d*NF+f} is the 
    observation index in {0..NP-1}, and {d} is the index of the signifcant 
    window position in {0..ND-1}. This is the dependent variable for the
    linear regression.

    Then  the   vector  {fr[0..NS-1]}   is  mapped  to   another  vector
    {gr[0..NS-1]}   by   the   orthonormal    linar   map   defined   by
    {multifok_focus_op_basis}. The  procedure computes a  certain number
    {NQ} of  quadratic terms {Q[0..NQ-1]} from  {gr[0..NS-1]}. It stores
    those terms in {(*termsP)[p*NQ + q]} for {q} in {[0..NQ-1]}.
    
    Furthermore, for each frame index {f} and each quadratic term index {q},
    the procedure also allocates a monochrome float image
    {qtimg[f*NQ + q]} and stores into it the value of {Q[q]} for each pixel of frame {f}.
    The vector {qtimg} of {NF*NQ} image pointers is also allocated by the
    procedure, and is returned in {*qtimgP}.
    
    The procedure also returns in {tnameP[0..NQ-1]} a newly allocated
    string that describes each term. In this description, "P{s}", for
    {s} in {0..NS-1}, denote the element {gr[s]} of the frame window
    vector after normalization and orthonormal mapping.
    
    The procedure also returns in {*ixP[d]} and {*iyP[d]} the 
    column and row of the center pixel of each significant window position.
    
    The procedure allocates the table {*tnameP} (with {NQ} elements),
    the matrix {*termsP} (with {NP} rows and {NQ} columns, linearized
    by rows), the vectors {*ixP,*iyP} (with {ND} elements),
    and the vectors {*e2P,*fminP} (with {NP} elements). It also
    returns {NP} in {*NPP}, {ND} in {*NDP}, and {NQ} in {NQP}.
    
    All averages, deviations, remappings, and squared distance
    computations use the window pixel weights defined by
    {multifok_focus_op_prod_weights}. The procedure assumes that all
    images are contaminated with Gaussian noise with zero mean and
    standard deviation {noise}. This mainly affects the normalization of
    the sample vectors. */
    
#endif
