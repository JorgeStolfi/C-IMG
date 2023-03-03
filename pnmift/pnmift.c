#define PROG_NAME "pnmift"
#define PROG_DESC "compute the Image Foresting Transform (IFT) of a PBM/PGM/PPM image"
#define PROG_VERS "1.0"

#define pnmift_C_COPYRIGHT "Copyright © 2001 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-02-12 07:53:24 by stolfi */

#define PROG_HELP \
  PROG_NAME  " \\\n" \
  "  -radius {RADIUS} \\\n" \
  "  -rootCost {ROOTFN} \\\n" \
  "  -arcCost {ARCFN} \\\n" \
  "  -pathCost {PATHFN} \\\n" \
  "  [ -seedImage {PGMFILE} ] [ -seedPixels {TXTFILE} ] \\\n" \
  "  { -output  {PREFIX} | \\\n" \
  "    {  -costs | -preds | -labels | -extract {LABEL} } {OUTPGM} \\\n" \
  "  } \\\n" \
  "  [ -roots  {OUTPGM} ] [ -spread {OUTPNM} ] \\\n" \
  "  [ -boxes {OUTTXT} [ -margin {MRG} ]] \\\n" \
  "  [ -plot {OUTEPS} ]  \\\n" \
  "  [ -bgColor {VALUE} {VALUE} {VALUE} | -bgGray {VALUE} ] \\\n" \
  "  [ -maxCostValue {MAXCOST} ] [ -maxCostPixel {MAXPIXEL} ] \\\n" \
  "  [ -lifo | -fifo ]  [ -reverse ] \\\n" \
  "  [ -verbose ] [ -listFunctions ] \\\n" \
  "  [ {PNMFILE} ]"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPT1 "\n" \
  "\n" \
  PROG_INFO_OPT2 "\n" \
  "\n" \
  PROG_INFO_OPT3 "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pgm(5), ppm(5), pnm(5), pgmselect(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in sep/2001 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "  IFT concept and algorithm by A.X.Falcão, R.Lotufo and J.Stolfi." \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2010-06-03 J. Stolfi:\n" \
  "    Substantial rewrite, using {libift} and generic" \
  " {pqueue.h} instead of {ift_queue.h}.  Costs may be arbitrary {double}s, but" \
  " time is now {N*log(N)} and partially broke \"-lifo\"/\"-fifo\".\n" \
  "\n" \
  "BUGS\n" \
  "  The options \"-lifo\"/\"-fifo\" do not affect queue order" \
  " (but still affect other tie-breaking).  Better tie-breaking" \
  " tools are needed.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmift_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program Reads a portable anymap {PNMFILE} and computes" \
  " its Image Foresting Transform (IFT), for a given set of seed" \
  " pixels and a specified path cost function.\n" \
  "\n" \
  "  The seed pixels can be specified in two ways: as a greyscale" \
  " image, or through a text file.  If neither method is used," \
  " every pixel is assumed to be a seed.\n" \
  "\n" \
  "  The IFT is written as three separate graymaps," \
  " containing the label, predecessor, and cost of each pixel," \
  " respectively.  Each image can be requested and named separately" \
  " with the \"-cost\", \"-pred\", and \"-label\" directives.  Alternatively," \
  " one may specify a common filename prefix as the \"-output {PREFIX}\" parameter;" \
  " the three images will be written out as \"{PREFIX}-costs.pgm\"," \
  " \"{PREFIX}-preds.pgm\", and \"{PREFIX}-labels.pgm\".\n" \
  "\n" \
  "  Optionally, the program may also writes a binary image showing" \
  " the root pixels (see the \"-roots\" directive); and/or a copy" \
  " of the original image where each root pixel has been replicated" \
  " throughout its influence zone (see the \"-spread\" directive); and/or" \
  " an image showing only the pixels with a specified label (see the \"-extract\" directive); and/or" \
  " a text file with the bounding box of each region (see the \"-boxes\" directive)."
  
#define PROG_INFO_OPT1 \
  "  -seedImage {PGMFILE}\n" \
  "    where {PGMFILE} is the name of a a grayscale image file in" \
  " the \"pgm\" format, where the non-zero pixels are the seeds," \
  " and their values are the respective labels.  The root cost of each pixel" \
  " is set to {+oo} if the label is zero, otherwise to the value specified" \
  " by the \"-rootCost\" argument.\n" \
  "\n" \
  "  -seedPixels {TEXTFILE}\n" \
  "    Specifies the name of an ascii file containing the seed pixels, their" \
  " labels, and (optionally) their root costs. Each line" \
  " of the file should have the format\n" \
  "      \"(\" {COL} {ROW} \")\" \"=\" {LABEL} [ {ROOTCOST} ]\n" \
  "    The {COL} and {ROW} are counted from the upper left corner, starting at 0.  The" \
  " label should be an integer in the range [0..2^16-1], and the {ROOTCOST}, if present, should" \
  " be a positive real number or \"INF\".  If {ROOTCOST} is omitted, it defaults to {+oo}" \
  " if the label is zero, otherwise to the value specified by the \"-rootCost\" argument.\n" \
  "    If both \"-seedImage\" and" \
  " \"-seedPixels\" are specified, the {PGMFILE} is read and used first, then" \
  " the {TEXTFILE} is read and its settings override those of the image file. If neither \"-seedImage\"" \
  " nor \"-seedPixels\"is specified, every pixel is" \
  " assigned label 1 and its root cost is set as specified by the \"-rootCost\" argument.\n" \
  "\n" \
  "  -radius {NUM}\n" \
  "    Specifies the radius of the neighborhood to use in the" \
  " propagation algorithm.  Use 1.0 for 4-neighbor topology, 1.5 for" \
  " 8-neighbor.  This parameter is required, and must be within 1.0" \
  " and 127.0, inclusive."
  
#define PROG_INFO_OPT2 \
  "  -pathCost {PATHFN}\n" \
  "    The name of the path-cost function to use.  The valid names are:\n" \
  pnmift_path_cost_fn_INFO "\n" \
  "\n" \
  "  -arcCost {ARCFN}\n" \
  "    The name of the arc-cost function to use.  The valid names are:\n" \
  pnmift_arc_cost_fn_INFO "\n" \
  "\n" \
  "  -output {PREFIX}\n" \
  "    The common prefix for output file names.  Use either this" \
  " option, or the options \"-costs\", \"-labels\", \"-preds\".\n" \
  "\n" \
  "  -costs {PGMFILE}\n" \
  "    Specifies the name of the output image file that contains" \
  " the pixel costs.\n" \
  "\n" \
  "  -preds {PGMFILE}\n" \
  "    Specifies the name of the output image file that contains" \
  " the predecessors of the pixels in the minimum-path forest.  Each" \
  " predecessor is given as two relative displacements (dh,dv), each" \
  " a signed integer in [-127..+127], packed into a 16-bit pixel" \
  " value 256*(dv+128) + (dh+128). \n" \
  "\n" \
  "  -labels {PGMFILE}\n" \
  "    Specifies the name of the output image file that contains" \
  " the propagated pixel labels.  Each label is stored as" \
  " a 16-bit pixel value.  A zero label marks those pixels" \
  " that cannot be reached from any seed at finite path cost.\n" \
  "\n" \
  "  -roots {PGMFILE}\n" \
  "    Causes the program to write the named binary output" \
  " image image (a graymap with maxval=1), showing which pixels" \
  " are roots (value 1) or non-roots (value 0).\n" \
  "\n" \
  "  -spread {PNMFILE}\n" \
  "    Causes the program to write an image of the same type and" \
  " depth as the original, where every pixel has been set to the" \
  " color of its correponding root pixel. \n" \
  "\n" \
  "  -extract {LABEL} {PNMFILE}\n" \
  "    Causes the program to write an image of the same type and" \
  " depth as the original, where every pixel of class {LABEL} is copied" \
  " from the original image, while all other pixels are set to the" \
  " background pixel value. \n" \
  "\n" \
  "  -boxes {TXTFILE}\n" \
  "    Causes the program to write a text file containg one" \
  " line for each possible label, with the bounding box" \
  " of the corresponding region, in the format\n" \
  "\n" \
  "      {LABEL} {HMIN} {VMIN} {HSIZE} {VSIZE}\n" \
  "\n" \
  "    If the region is empty, the last four fields are zero."
  
#define PROG_INFO_OPT3 \
  "  -margin {MARGIN}\n" \
  "    Used in conjuntion with \"-boxes\".  Expands the bounding" \
  " box of each non-empty region by the given number of pixels, in" \
  " all four directions, then clips the result against the image bounds.\n" \
  "\n" \
  "  -bgColor {VALUE} {VALUE} {VALUE}\n" \
  "  -bgGray {VALUE}\n" \
  "    These options specify the background pixel value to use for" \
  " the \"-extract\" option.  Only one of the two should be" \
  " specified.  If \"-bgColor\" is used for gray" \
  " image input, only the first {VALUE} is used.  The default" \
  " is \"-bgColor 0 0 0\" = \"-bgGray 0\".\n" \
  "\n" \
  "  -plot {EPSFILE}\n" \
  "    Causes the program to write an Encapsulated Postscript version" \
  " of the image, overlaid with a drawing of the minimum-path forest.\n" \
  "\n" \
  "  -maxCostPixel {MAXPIXEL}\n" \
  "    Specifies the maximum pixel value to be used when writing the" \
  " cost layer of the IFT.  Defaults to 2^16-1 = 65535.\n" \
  "\n" \
  "  -maxCostValue {MAXCOST}\n" \
  "    Specifies the maximum nominal path cost, for the purpose of" \
  " output scaling.  The units are the same used internally by the" \
  " IFT procedure (see {ift.h}). Defaults to the maximum path cost" \
  " observed in the IFT.  When writing the cost layer, infinte path" \
  " costs will be mapped to MAXPIXEL.  Finite path costs in" \
  " [0..MAXCOST] will be mapped linearly to [0..{MAXPIXEL-1} ].   Finite" \
  " costs higher than {MAXCOST} will be mapped to {MAXPIXEL-1}.\n" \
  "\n" \
  "  -lifo\n" \
  "    Specifies LIFO precedence for selecting between equal-cost paths.\n" \
  "\n" \
  "  -fifo\n" \
  "    Specifies FIFO precedence for selecting between equal-cost paths. \n" \
  "\n" \
  "  -reverse\n" \
  "    Reverses the scanning order of seeds and outgoing arcs.\n" \
  "\n" \
  "  -verbose\n" \
  "    Produces diagnostic output."

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <fget.h>
#include <affirm.h>
#include <argparser.h>
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>
#include <jspnm.h>
#include <jsstring.h>
#include <jsfile.h>
#include <epswr.h>

#include <ift.h>
#include <ift_plot.h>
#include <pnmift_root_cost_fn.h>
#include <pnmift_arc_cost_fn.h>
#include <pnmift_path_cost_fn.h>

#include <pnmift_image.h>

/* PROTOTYPES */

typedef struct options_t
  { /* Input image: */
    char *input_img_name;           
    /* Seed sources: */
    char *seed_img_name; 
    char *seed_txt_name; 
    /* Adjacency relation: */
    double radius;
    /* Path costs: */
    pnmift_root_cost_fn_t *rfn;
    pnmift_arc_cost_fn_t *afn;
    pnmift_path_cost_fn_t *pfn;
    /* Tie breaking: */
    ift_tie_breaking_t tbreak;
    ift_scan_order_t order;
    /* Output images: */
    uint16_t bgColor[3];
    ift_path_cost_t max_cost_value;
    uint16_t max_cost_pixel;
    char *cost_img_name;
    char *pred_img_name;
    char *label_img_name;
    char *root_img_name;
    char *spread_img_name;
    char *extract_img_name;
    uint16_t extract_label;
    char *region_boxes_name;
    int32_t box_margin;
    char *plot_eps_name; 
    /* Diagnostics: */
    int32_t verbose;
  } options_t;
  /* Command line arguments. */

int32_t main(int32_t argc, char* argv[]);

options_t *pnmift_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments. */

void pnmift_read_image_as_graph(options_t *o, uint16_t *input_maxvalP, int32_t *input_chnsP, ift_graph_t **GP, frgb_t **rgbP);
  /* Reads the PNM file "{o->input_img_name}" and transforms it into
    an IFT graph, stored in {*GP}.  Returns the image pixels in {*rgbP}, indexed
    by the graphnode indices.  Also updates {*input_maxvalP} to be at least as big
    as the nominal maximum value of any sample in the image (the PNM
    {maxval} parameter). */
    
void pnmift_set_all_labels(ift_graph_t *G, uint16_t label[], uint16_t lab);
  /* Sets {label[k]} to {lab} for all {k} in {0..G->nodes-1}. */

void pnmift_set_root_costs_from_labels
  ( ift_graph_t *G, 
    frgb_t rgb[], 
    int32_t chns, 
    uint16_t label[],
    pnmift_root_cost_fn_t *rfn
  );
  /* Computes the root costs (trivial path costs, handicaps) {G->node[i].C}
    for all nodes.  If {label[i]} is {pnmift_NON_SEED_LABEL},
    sets the cost to {+oo}, otherwise sets it to {rfn(rgb[i],chns)}. */

void pnmift_set_labels_from_image_file
  ( options_t *o, 
    ift_graph_t *G, 
    uint16_t label[]
  );
  /* Sets the labels from the PNM image file named
    "{o->seed_img_name}". Also updates {*maxvalP} to to be at least as
    big as the nominal maximum value of any sample in the image (the
    PNM {maxval} parameter). */

void pnmift_set_labels_and_costs_from_text_file
  ( options_t *o, 
    ift_graph_t *G, 
    frgb_t rgb[], 
    int32_t chns, 
    uint16_t label[],
    pnmift_root_cost_fn_t *rfn
  );
  /* Sets the labels {label[0..G.nodes-1]} and possibly the root costs
    {G->node[0..G.nodes-1].C} from the text file named
    "{o->seed_txt_name}". */

void pnmift_compute_forest
  ( ift_graph_t *G, 
    frgb_t rgb[],
    int32_t chns,
    pnmift_arc_cost_fn_t *afn, 
    pnmift_path_cost_fn_t *pfn, 
    ift_tie_breaking_t tbreak, 
    ift_scan_order_t order, 
    ift_path_cost_t *maxCostP
  );
  /* Computes the IFT on the image graph {G}. Assumes that all paths
    in {G} are trivial and that the costs of those trivial paths (the
    root costs) have been stroed in {G.node[0..G.nodes-1].C}. The cost
    of a path is derived from those root costs and the arc costs,
    combined recursively along the path by the function {pfn}. The
    cost of an arc {(p,q)} is presumably computed by {afn} from the
    image values {rgb[p],rgb[q]} and the arc displacement {q - p}. */

void pnmift_write_cost_image(options_t *o, ift_graph_t *G, ift_path_cost_t actual_max_cost);
void pnmift_write_label_image(options_t *o, ift_graph_t *G, uint16_t label[]);
void pnmift_write_pred_image(options_t *o, ift_graph_t *G);
void pnmift_write_root_image(options_t *o, ift_graph_t *G);
void pnmift_write_spread_image(options_t *o, ift_graph_t *G, frgb_t rgb[], int32_t chns, uint16_t maxval);
void pnmift_write_extract_image(options_t *o, ift_graph_t *G, uint16_t label[], frgb_t rgb[], int32_t chns, uint16_t maxval);
  /* These procedure construct various output images (by the corresponding procedures
     in {ift.h}) and write them to files called "{o->cost_img_name}",
     "{label_img_name}", etc.. */

void pnmift_write_region_bboxes(options_t *o, ift_graph_t *G, uint16_t label[]);
  /* Writes a text file containing a rectangular bounding box for each 
    region of {G} in {0..maxlabel}, where {maxlabel} is the maximum label
    appearing in {label[0..G.nodes]}.  See {ift_write_boxes} in {ift.h}. */

void pnmift_write_postscript_plot(options_t *o, ift_graph_t *G, frgb_t rgb[]);
  /* Writes an Encapsulated Postscript file with a plot of the 
    spanning forest contained in {G}. */

void skip_rest_of_line(FILE *f);
  /* Skip blanks/comments until end of line, complain if anything else: */

/* REAL CODE */

int32_t main(int32_t argc, char** argv)
  {
    /* Parse command line arguments. */
    options_t *o = pnmift_parse_options(argc, argv);

    int32_t input_chns;                  /* Channels in input image (1 or 3). */
    uint16_t input_maxval;       /* Max sample value of input image. */
    ift_path_cost_t actual_max_cost; /* Actual maximum cost observed in the IFT. */

    /* Get input image, and turn it into a graph {G}: */
    ift_graph_t *G;
    frgb_t *rgb;
    pnmift_read_image_as_graph(o, &input_maxval, &input_chns, &G, &rgb); /* Graph of input image. */
    
    /* Set labels and root costs for seed pixels: */
    uint16_t *label = notnull(malloc(G->nodes*sizeof(uint16_t)), "out of memory");
    if ((o->seed_img_name == NULL) && (o->seed_txt_name == NULL)) 
      { /* No seed pixels are specified; mark all pixels as seeds: */
        if (o->verbose) { fprintf(stderr, "setting all seed labels to %d...\n", pnmift_DEFAULT_SEED_LABEL); }
        pnmift_set_all_labels(G, label, pnmift_DEFAULT_SEED_LABEL);
        /* Set root costs of seed pixels: */
        pnmift_set_root_costs_from_labels(G, rgb, input_chns, label, o->rfn);
      }
    else
      { /* Mark the seed pixels given as a PNM image: */
        if (o->seed_img_name != NULL)
          { pnmift_set_labels_from_image_file(o, G, label); }
        else
          { pnmift_set_all_labels(G, label, pnmift_NON_SEED_LABEL); }
        pnmift_set_root_costs_from_labels(G, rgb, input_chns, label, o->rfn);
        /* Modify the seed pixels given as a text file: */
        if (o->seed_txt_name != NULL)
          { if (o->verbose) { fprintf(stderr, "reading seed file %s...\n", o->seed_txt_name); }
            pnmift_set_labels_and_costs_from_text_file(o, G, rgb, input_chns, label, o->rfn);
          }
      }

    /* Compute the IFT: */
    if (o->verbose) { fprintf(stderr, "computing the IFT...\n"); }
    pnmift_compute_forest(G, rgb, input_chns, o->afn, o->pfn, o->tbreak, o->order, &actual_max_cost);
      
    /* Write the output images: */
    
    /* Cost layer: */
    if (o->cost_img_name != NULL)
      { pnmift_write_cost_image(o, G, actual_max_cost); }
      
    /* Label layer: */
    if (o->label_img_name != NULL)
      { pnmift_write_label_image(o, G, label); }
    
    /* Predecessor layer: */
    if (o->pred_img_name != NULL)
      { pnmift_write_pred_image(o, G); }
      
    /* Roots indicator file: */
    if (o->root_img_name != NULL)
      { pnmift_write_root_image(o, G); }
      
    /* Root-spread image: */
    if (o->spread_img_name != NULL)
      { pnmift_write_spread_image(o, G, rgb, input_chns, input_maxval); }
      
    /* Single-label image: */
    if (o->extract_img_name != NULL)
      { pnmift_write_extract_image(o, G, label, rgb, input_chns, input_maxval); }

    /* Region bounding boxes: */
    if (o->region_boxes_name != NULL)
      { pnmift_write_region_bboxes(o, G, label); }
      
    /* Postscript plot: */
    if (o->plot_eps_name != NULL)
      { pnmift_write_postscript_plot(o, G, rgb); }
      
    /* All done: */
    if (o->verbose) { fprintf(stderr, "done.\n"); }
    exit(0);
  }

void pnmift_compute_forest
  ( ift_graph_t *G, 
    frgb_t rgb[],
    int32_t chns,
    pnmift_arc_cost_fn_t *afn, 
    pnmift_path_cost_fn_t *pfn, 
    ift_tie_breaking_t tbreak, 
    ift_scan_order_t order, 
    ift_path_cost_t *maxCostP
  )
  {
    auto ift_path_cost_t path_cost(ift_node_t *s, ift_rel_arc_t *a, ift_node_t *t, ift_graph_t *G);
    
    ift_compute_forest(G, path_cost, tbreak, order, maxCostP);
    return;
    
    ift_path_cost_t path_cost(ift_node_t *s, ift_rel_arc_t *a, ift_node_t *t, ift_graph_t *G)
      {
        demand((s == NULL) == (a == NULL), "s/a inconstency");
        demand(t != NULL, "t is null");
        if (s == NULL)
          { return t->C; }
        else
          { ift_node_index_t si = ift_node_index(G, s->col, s->row);
            ift_node_index_t ti = ift_node_index(G, t->col, t->row);
            double aC = afn(rgb[si], a, rgb[ti], chns);
            return pfn(s->C, aC);
          }
      }
  }

void pnmift_read_image_as_graph
  ( options_t *o, 
    uint16_t *input_maxvalP, 
    int32_t *input_chnsP, 
    ift_graph_t **GP, 
    frgb_t **rgbP
  )
  { uint16_image_t *input_img;
    input_img = uint16_image_read_pnm_named(o->input_img_name, o->verbose);
    if (o->verbose) { fprintf(stderr, "maxval = %d\n", input_img->maxval); }
    if (o->verbose) { fprintf(stderr, "building image graph...\n"); }
    ift_graph_t *G = ift_make_graph(input_img->cols, input_img->rows, o->radius);
    frgb_t *rgb = notnull(malloc(G->nodes*sizeof(frgb_t)), "no mem");
    if (o->verbose) { fprintf(stderr, "cols = %d rows = %d\n", G->cols, G->rows); }
    if (o->verbose) { fprintf(stderr, "channels = %d\n", input_img->chns); }
    if (o->verbose) { fprintf(stderr, "nodes = %d arcs/node = %d\n", G->nodes, G->arcs); }
    pnmift_set_values_from_image(input_img, G, rgb);
    (*input_maxvalP) = input_img->maxval;
    (*input_chnsP) = input_img->chns;
    (*GP) = G;
    (*rgbP) = rgb;
    uint16_image_free(input_img);
  }
    
void pnmift_set_all_labels(ift_graph_t *G, uint16_t label[], uint16_t lab)
  { uint64_t LL = lab; /* To supress a bogus "always false" gcc warning. */
    demand(LL <= pnmift_MAX_LABEL, "bad label");
    int32_t i;
    for (i = 0; i < G->nodes; i++) { label[i] = lab; }
  }
  
void pnmift_set_root_costs_from_labels
  ( ift_graph_t *G, 
    frgb_t rgb[], 
    int32_t chns, 
    uint16_t label[],
    pnmift_root_cost_fn_t *rfn
  )
  { int32_t i;
    for (i = 0; i < G->nodes; i++) 
      { if (label[i] == pnmift_NON_SEED_LABEL)
          { G->node[i].C = +INF; }
        else
          { G->node[i].C = rfn(rgb[i], chns); }
      }
  }

void pnmift_set_labels_from_image_file
  ( options_t *o, 
    ift_graph_t *G, 
    uint16_t label[]
  )
  { uint16_image_t *seed_img;
    seed_img = uint16_image_read_pnm_named(o->seed_img_name, o->verbose);
    if (seed_img->chns != 1) { pnm_error("seed image must be monochromatic"); }
    if (o->verbose) { fprintf(stderr, "marking seeds in graph...\n"); }
    pnmift_set_labels_from_image(seed_img, G, label, o->verbose);
    uint16_image_free(seed_img);
  }

void pnmift_set_labels_and_costs_from_text_file
  ( options_t *o, 
    ift_graph_t *G, 
    frgb_t rgb[], 
    int32_t chns, 
    uint16_t label[],
    pnmift_root_cost_fn_t *rfn
  )
  {
    FILE *f = open_read(o->seed_txt_name, o->verbose);
    int32_t ch;
    ch = fgetc(f);
    while (ch != EOF)
      { if ((ch == ' ') || (ch == '\n') || (ch == '\011') || (ch == '\015'))
          { /* Ignore. */ }
        else if (ch == '#')
          { while ((ch != '\n') && (ch != EOF)) { ch = fgetc(f); } }
        else if (ch == '(')
          { int32_t col = fget_int32(f); 
            demand((col >= 0) && (col < G->cols), "bad col in seed file"); 
            int32_t row = fget_int32(f);
            demand((row >= 0) && (row < G->rows), "bad row in seed file");
            fget_skip_spaces_and_match(f, ")");
            fget_skip_spaces_and_match(f, "=");
            int32_t lab = fget_int32(f);
            demand((lab >= 0) && (lab <= pnmift_MAX_LABEL), "bad label in seed file");
            int32_t ig = ift_node_index(G, (int16_t)col, (int16_t)row);
            ift_node_t *pg = &(G->node[ig]);
            assert(pg->col == col);
            assert(pg->row == row);
            label[ig] = (uint16_t)lab;
            /* Read handicap if any, else uses handicap = 0: */
            fget_skip_spaces(f);
            int32_t ch1 = fgetc(f);
            if (ch1 != EOF) { ungetc(ch1, f); }
            if ((ch1 == '+') || (ch1 == '-') || (ch1 == '.') || ((ch1 >= '0') && (ch1 <= '9')))
              { pg->C = fget_double(f); }
            else
              { pg->C = rfn(rgb[ig], chns); }
            skip_rest_of_line(f);
            if (o->verbose)
              { fprintf(stderr, "  ( %5d %5d ) = %d %f\n", col, row, label[ig], pg->C); }
          }
        else
          { fprintf(stderr, "ch = «%c»\n", ch);
            demand(FALSE, "unexpected char in seed file");
          }
        if (ch != EOF) { ch = fgetc(f); }
      }
    fclose(f);
  }

void pnmift_write_cost_image(options_t *o, ift_graph_t *G, ift_path_cost_t actual_max_cost)
  { if (o->verbose) { fprintf(stderr, "extracting the cost image...\n"); }
    fprintf(stderr, "maximum observed path cost = %12.8f\n",  actual_max_cost);
    if (o->max_cost_value == 0) { o->max_cost_value = actual_max_cost; }
    if (o->verbose) 
      { fprintf
         ( stderr, "scaling [0..%.0f] to [0..%d], +oo to %d\n",
           o->max_cost_value, 
           o->max_cost_pixel-1,
           o->max_cost_pixel
         );
      }
    uint16_image_t *cost_img = pnmift_get_cost_image(G, o->max_cost_value, o->max_cost_pixel);
    uint16_image_write_pnm_named(o->cost_img_name, cost_img, 0, o->verbose);
    uint16_image_free(cost_img);
  }

void pnmift_write_label_image(options_t *o, ift_graph_t *G, uint16_t label[])
  { if (o->verbose) { fprintf(stderr, "extracting the label image...\n"); }
    /* Find maximum label (minimum 1): */
    uint16_t maxlabel = 1;
    int32_t i;
    for (i = 0; i < G->nodes; i++) { if (label[i] > maxlabel) { maxlabel = label[i]; } }
    uint16_image_t *label_img = pnmift_get_label_image(G, label, maxlabel);
    uint16_image_write_pnm_named(o->label_img_name, label_img, 0, o->verbose);
    uint16_image_free(label_img);
  }

void pnmift_write_pred_image(options_t *o, ift_graph_t *G)
  { 
    if ((int32_t)PNM_FILE_MAX_MAXVAL < 65535) 
      { fprintf(stderr, "PNM_FILE_MAX_MAXVAL = %d too small - pred image not generated\n", PNM_FILE_MAX_MAXVAL); }
    else
      { if (o->verbose) { fprintf(stderr, "extracting the pred image...\n"); }
        uint16_image_t *pred_img = pnmift_get_pred_image(G);
        uint16_image_write_pnm_named(o->pred_img_name, pred_img, 0, o->verbose);
        uint16_image_free(pred_img);
      }
  }

void pnmift_write_root_image(options_t *o, ift_graph_t *G)
  { if (o->verbose) { fprintf(stderr, "generating the root mask...\n"); }
    uint16_image_t *root_img = pnmift_get_root_image(G);
    uint16_image_write_pnm_named(o->root_img_name, root_img, 0, o->verbose);
    uint16_image_free(root_img);
  }

void pnmift_write_spread_image(options_t *o, ift_graph_t *G, frgb_t rgb[], int32_t chns, uint16_t maxval)
  { uint16_image_t *spread_img;
    if (o->verbose) { fprintf(stderr, "extracting the spread image...\n"); }
    spread_img = pnmift_get_spread_image(G, rgb, chns, maxval);
    uint16_image_write_pnm_named(o->spread_img_name, spread_img, 0, o->verbose);
    uint16_image_free(spread_img);
  }

void pnmift_write_extract_image(options_t *o, ift_graph_t *G, uint16_t label[], frgb_t rgb[], int32_t chns, uint16_t maxval)
  { if (o->verbose) { fprintf(stderr, "extracting the sub-image with label = %d...\n", o->extract_label); }
    uint16_image_t *extract_img = pnmift_get_single_label_image(G, label, o->extract_label, rgb, o->bgColor, chns, maxval);
    if (o->verbose) { fprintf(stderr, "maxval = %d\n", extract_img->maxval); }
    uint16_image_write_pnm_named(o->extract_img_name, extract_img, 0, o->verbose);
    uint16_image_free(extract_img);
  }

void pnmift_write_region_bboxes(options_t *o, ift_graph_t *G, uint16_t label[])
  { FILE *wr = open_write(o->region_boxes_name, o->verbose);
    if (o->verbose) { fprintf(stderr, "writing the bounding boxes...\n"); }
    /* Find maximum seed label: */
    uint16_t maxlabel = 0;
    int32_t i;
    for (i = 0; i < G->nodes; i++) { if (label[i] > maxlabel) { maxlabel= label[i]; } }
    pnmift_write_boxes(wr, G, label, maxlabel, o->box_margin);
    fclose(wr);
  }

void pnmift_write_postscript_plot(options_t *o, ift_graph_t *G, frgb_t rgb[])
  { double hsize = 12.0*(double)G->cols;
    double vsize = 12.0*(double)G->rows;
    if (o->verbose) { fprintf(stderr, "generating the tree plot %s...\n", o->plot_eps_name); }
    FILE *psfile = open_write(o->plot_eps_name, o->verbose);
    epswr_figure_t *eps = epswr_new_figure(psfile, hsize, vsize, 4.0, 4.0, 4.0, 4.0, FALSE);
    epswr_set_window
      ( eps,
        2.0, 2.0 + hsize, 2.0, 2.0 + vsize, 
        FALSE,
        0.0, (double)G->cols, 0.0, (double)G->rows
      );
    /* epswr_set_grid(eps, G->cols, G->rows); */
    ift_plot_pixel_values(eps, G, rgb, 0.5);
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.5, 0.0, 0.0); 
    ift_plot_forest_edges(eps, G);
    epswr_set_pen(eps, 0.0,0.0,0.0, 0.25, 0.0, 0.0); 
    frgb_t root_rgb = (frgb_t){{ 1,1,1 }};
    ift_plot_forest_nodes(eps, G, TRUE,  1.0,  &root_rgb);
    frgb_t plain_rgb = (frgb_t){{ 0,0,0 }};
    ift_plot_forest_nodes(eps, G, FALSE, 0.5,  &plain_rgb);
    epswr_end_figure(eps);
  }

options_t *pnmift_parse_options(int32_t argc, char **argv)
  {
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    if (argparser_keyword_present(pp, "-seedImage"))
      { o->seed_img_name = argparser_get_next(pp); }
    else 
      { o->seed_img_name = NULL; }
    
    if (argparser_keyword_present(pp, "-seedPixels"))
      { o->seed_txt_name = argparser_get_next(pp); }
    else  
      { o->seed_txt_name = NULL; }
    
    if (argparser_keyword_present(pp, "-radius"))
      { o->radius = argparser_get_next_double(pp, 1.0, 127.0); }
    else  
      { argparser_error(pp, "must specify a radius"); }
    
    if (argparser_keyword_present(pp, "-rootCost"))
      { char *rfn_str = argparser_get_next(pp);
        o->rfn = pnmift_root_cost_fn_from_name(rfn_str);
      }
    else  
      { argparser_error(pp, "must specify a root cost function"); }
    
    if (argparser_keyword_present(pp, "-arcCost"))
      { char *afn_str = argparser_get_next(pp);
        o->afn = pnmift_arc_cost_fn_from_name(afn_str);
      }
    else  
      { argparser_error(pp, "must specify an arc cost function"); }
    
    if (argparser_keyword_present(pp, "-pathCost"))
      { char *pfn_str = argparser_get_next(pp);
        o->pfn = pnmift_path_cost_fn_from_name(pfn_str);
      }
    else  
      { argparser_error(pp, "must specify a path cost function"); }
    
    if (argparser_keyword_present(pp, "-output"))
      { /* Define matched filenames for costs, labels, and preds images: */
        char *output_str = argparser_get_next(pp);
        o->cost_img_name = txtcat(output_str, "-costs.pgm"); 
        o->label_img_name = txtcat(output_str, "-labels.pgm"); 
        o->pred_img_name = txtcat(output_str, "-preds.pgm"); 
      }
    else  
      { /* Allow individual filenames for costs, labels, and preds images: */
        if (argparser_keyword_present(pp, "-costs"))
          { o->cost_img_name = argparser_get_next(pp); }
        else  
          { o->cost_img_name = NULL; }

        if (argparser_keyword_present(pp, "-preds"))
          { o->pred_img_name = argparser_get_next(pp); }
        else  
          { o->pred_img_name = NULL; }

        if (argparser_keyword_present(pp, "-labels"))
          { o->label_img_name = argparser_get_next(pp); }
        else  
          { o->label_img_name = NULL; }
      }
    
    if (argparser_keyword_present(pp, "-roots"))
      { o->root_img_name = argparser_get_next(pp); }
    else  
      { o->root_img_name = NULL; }
    
    if (argparser_keyword_present(pp, "-spread"))
      { o->spread_img_name = argparser_get_next(pp); }
    else  
      { o->spread_img_name = NULL; }
    
    if (argparser_keyword_present(pp, "-extract"))
      { o->extract_label = (uint16_t)argparser_get_next_int(pp, 0, INT32_MAX);
        o->extract_img_name = argparser_get_next(pp);
      }
    else
      { o->extract_label = 0;
        o->extract_img_name = NULL; 
      }

    if (argparser_keyword_present(pp, "-boxes"))
      { o->region_boxes_name = argparser_get_next(pp); }
    else  
      { o->region_boxes_name = NULL; }
      
    if 
      ( (o->cost_img_name == NULL) &&
        (o->pred_img_name == NULL) &&
        (o->label_img_name == NULL) &&
        (o->root_img_name == NULL) &&
        (o->spread_img_name == NULL) &&
        (o->region_boxes_name == NULL)
      )
      { fprintf(stderr, "what, no output files?\n"); }
    
    if (argparser_keyword_present(pp, "-margin"))
      { o->box_margin = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else  
      { o->box_margin = 0; }
    
    if (argparser_keyword_present(pp, "-plot"))
      { o->plot_eps_name = argparser_get_next(pp); }
    else  
      { o->plot_eps_name = NULL; }

    if (argparser_keyword_present(pp, "-bgColor"))
      { int32_t ch;
        for (ch = 0; ch < 3; ch++) 
          { o->bgColor[ch] = (uint16_t)argparser_get_next_int(pp, 0, INT32_MAX); }
      }
    else if (argparser_keyword_present(pp, "-bgGray"))
      { uint16_t val = (uint16_t)argparser_get_next_int(pp, 0, INT32_MAX);
        int32_t ch;
        for (ch = 0; ch < 3; ch++) { o->bgColor[ch] = val; }
      }
    else
      { int32_t ch;
        for (ch = 0; ch < 3; ch++) { o->bgColor[ch] = 0; }
      }

    if (argparser_keyword_present(pp, "-maxCostValue"))
      { o->max_cost_value = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX); }
    else  
      { o->max_cost_value = 0; /* 0 means actual max */ }
    
    if (argparser_keyword_present(pp, "-maxCostPixel"))
      { o->max_cost_pixel = (uint16_t)argparser_get_next_int(pp, 1, PNM_FILE_MAX_MAXVAL); }
    else  
      { o->max_cost_pixel = PNM_FILE_MAX_MAXVAL; }
    
    if (argparser_keyword_present(pp, "-lifo"))
      { o->tbreak = ift_tbreak_LIFO; }
    else if (argparser_keyword_present(pp, "-fifo"))
      { o->tbreak = ift_tbreak_FIFO; }
    else 
      { o->tbreak = ift_tbreak_FIFO; }
    
    if (argparser_keyword_present(pp, "-reverse"))
      { o->order = ift_order_DN; }
    else  
      { o->order = ift_order_UP; }
    
    o->verbose = argparser_keyword_present(pp, "-verbose");
    
    /* Get optional input file name: */
    argparser_skip_parsed(pp);
    if (argparser_next(pp) != NULL) 
      { o->input_img_name = argparser_get_next(pp); }
    else
      { o->input_img_name = "-"; }

    /* Check for extraneous arguments: */
    argparser_finish(pp);
    
    return o;
  }

void skip_rest_of_line(FILE *f)
  {
    fget_skip_spaces(f);
    if (fget_test_char(f, '#'))
      { int32_t r; char c;
        do { r = fgetc(f); c = (char)r; } while ((r != EOF) || (c != '\n')); }
    else
      { fget_eol(f); }
  }
