#define PROG_NAME "pnmclassif"
#define PROG_DESC "fuzzy segmentation of a PBM/PGM/PPM image by NN pixel classification"
#define PROG_VERS "1.0"

#define pnmclassif_C_COPYRIGHT "Copyright © 2010 by the State University of Campinas (UNICAMP)"

/* Last edited on 2023-10-14 23:06:15 by stolfi */

#define PROG_HELP \
  PROG_NAME  " \\\n" \
  "  -distance {ARCFN} \\\n" \
  "  { -seedImage {PGMFILE} | -seedPixels {TXTFILE} } \\\n" \
  "  -output  {PREFIX} \\\n" \
  "  [ -moments {OUTTXT} ] \\\n" \
  "  [ -verbose ] [ -listDistances ] \\\n" \
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
  "  Created in jun/2010 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "  IFT concept and algorithm by A.X.Falcão, R.Lotufo and J.Stolfi." \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  Option bla bla added by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmclassif_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program reads a portable anymap {PNMFILE} and outputs a set of grayscale masks" \
  " based on a local probabilistic classification of pixel colors into" \
  " a specified number {N} of classes, where each" \
  " class is determined by a user-given subset of `representative pixels' (`seeds').\n" \
  "\n" \
  "  The seed pixels can be specified in two ways: as a greyscale" \
  " image, or through a text file.\n" \
  "\n" \
  "  The output is a set of {N} separate graymaps," \
  " containing the probability that each pixel belongs to a specific" \
  " class.  The image files are named \"{PREFIX}-{K}.pgm\", where {K} is" \
  " a three-digit class index from 1 to {N}."
  
#define PROG_INFO_OPT1 \
  "  -seedImage {PGMFILE}\n" \
  "    where {PGMFILE} is the name of a a grayscale image file in" \
  " the \"pgm\" format, where the non-zero pixels are the seeds," \
  " and their values in {PGMFILE} are the respective classes.  Every seed" \
  " will have an assigned weight of 1.\n" \
  "\n" \
  "  -seedPixels {TEXTFILE}\n" \
  "    Specifies the name of an ascii file containing the seed pixels, their" \
  " classes, and (optionally) their root costs. Each line" \
  " of the file should have the format\n" \
  "      \"(\" {COL} {ROW} \")\" \"=\" {CLASS} [ {WEIGHT} ]\n" \
  "    The {COL} and {ROW} are counted from the upper left corner, starting at 0.  The" \
  " class should be an integer in the range [0..2^16-1], and the {WEIGHT}, if present, should" \
  " be a finite non-negative real number; if omitted, it defaults to {1}.  If {CLASS} is" \
  " positive, it means that pixel {COL,ROW} is a representative of the class {CLASS}.  Entries" \
  " with {CLASS=0} or {WEIGHT=0} are ignored.\n" \
  "    If both \"-seedImage\" and" \
  " \"-seedPixels\" are specified, the {PGMFILE} is read and used first, then" \
  " the {TEXTFILE} is read and its settings override those of the image file. If neither \"-seedImage\"" \
  " nor \"-seedPixels\"is specified, every pixel is" \
  " assigned class 1 and its root cost is set as specified by the \"-rootCost\" argument."
  
#define PROG_INFO_OPT2 \
  "  -distance {DISTFN}\n" \
  "    The name of the color distance function to use.  The valid names are:\n" \
  pnmclassif_dist_fn_INFO "\n" \
  "\n" \
  "  -output {PREFIX}\n" \
  "    The common prefix for output file names.\n" \
  "\n" \
  "  -moments {TXTFILE}\n" \
  "    Causes the program to write a text file containg one" \
  " line for each class, with the mean position, orientation and" \
  " principal moments of the class probabilities, in the format\n" \
  "\n" \
  "      {CLASS} {HMEAN} {VMEAN} {TILT} {USIGMA} {VSIGMA}\n" \
  "\n" \
  "    If the probability distribution is zero, the last five fields are zero."
  
#define PROG_INFO_OPT3 \
  "  -verbose\n" \
  "    Produces diagnostic output."

#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <fget.h>
#include <affirm.h>
#include <argparser.h>
#include <uint16_image.h>
#include <jspnm.h>
#include <jsstring.h>
#include <jsfile.h>
#include <pswr.h>

#include <pnmclassif_dist_fn.h>
#include <pnmclassif_image.h>

/* PROTOTYPES */

typedef struct options_t
  { /* Input image: */
    char *input_img_name;           
    /* Seed sources: */
    char *seed_img_name; 
    char *seed_txt_name; 
    /* Path costs: */
    pnmclassif_dist_fn_t *dfn;
    /* Output images: */
    char *dist_img_name;
    char *prob_img_name;
    char *class_moments_name;
    /* Diagnostics: */
    int verbose;
  } options_t;
  /* Command line arguments. */

int main(int argc, char* argv[]);

options_t *pnmclassif_parse_options(int argc, char **argv);
  /* Parses the command line arguments. */

float_image_t *pnmclassif_read_pnm_input_image(char *filename);
  /* Reads the PNM file "{filename}" and returns it as a {float_image_t}. */
    
void pnmclassif_get_seeds_from_image_file
  ( char *filename, 
    float_image_t *img, 
    int *nseedsP, 
    int *nclassesP,
    float **seedP, 
    double **weightP,
    int **classP
  );
  /* Read the PGM image file "{filename}".  For every pixel 
    from the PNM image file named
    "{o->seed_img_name}". Also updates {*maxvalP} to to be at least as
    big as the nominal maximum value of any sample in the image (the
    PNM {maxval} parameter). */

void pnmclassif_get_seeds_and_weights_from_text_file(char *filename, double **seedP, int *nseedsP);
  /* Sets the classes {class[0..G.nodes-1]} and possibly the root costs
    {G->node[0..G.nodes-1].C} from the text file named
    "{o->seed_txt_name}". */

float_image_t *pnmclassif_compute_classification
  ( float_image_t *img, 
    int nseeds, 
    int nclasses,
    float seed[], 
    double weight[],
    int class[]
  );
  /* Computes a pixel classification for the image {img}
    using seeds {seed[0..chns*nseeds-1]} with weights {0..nseeds-1]} and
    classes {class[0..nseeds-1]}, where {chns = img.sz[0]} is the number of 
    channels in the image {img}.  The classification 
    is returned as a float image with {nclasses+1} channels,
    where sample {x,y} of channel {c} is the probability that pixel
    {x,y} of {img} belongs to class {c}, for {c=1..nclasses}.
    !!! What about channel 0? !!! */

void pnmclassif_write_cost_image(options_t *o, ift_graph_t *G, ift_path_cost_t actual_max_cost);
void pnmclassif_write_class_image(options_t *o, ift_graph_t *G, uint16_t class[]);
void pnmclassif_write_pred_image(options_t *o, ift_graph_t *G);
void pnmclassif_write_root_image(options_t *o, ift_graph_t *G);
void pnmclassif_write_spread_image(options_t *o, ift_graph_t *G, frgb_t rgb[], int chns, uint16_t maxval);
void pnmclassif_write_extract_image(options_t *o, ift_graph_t *G, uint16_t class[], frgb_t rgb[], int chns, uint16_t maxval);
  /* These procedure construct various output images (by the corresponding procedures
     in {ift.h}) and write them to files called "{o->dist_img_name}",
     "{class_img_name}", etc.. */

void pnmclassif_write_class_moments(options_t *o, ift_graph_t *G, uint16_t class[]);
  /* Writes a text file containing a rectangular bounding box for each 
    region of {G} in {0..maxclass}, where {maxclass} is the maximum class
    appearing in {class[0..G.nodes]}.  See {ift_write_boxes} in {ift.h}. */

void skip_rest_of_line(FILE *f);
  /* Skip blanks/comments until end of line, complain if anything else: */

/* REAL CODE */

int main(int argc, char** argv)
  {
    /* Parse command line arguments. */
    options_t *o = pnmclassif_parse_options(argc, argv);

    int input_chns;                  /* Channels in input image (1 or 3). */
    uint16_t input_maxval;       /* Max sample value of input image. */

    /* Get input image: */
    float_image_t *img = pnmclassif_read_pnm_input_image(...);
    
    /* get seed pixels: */
    uint16_t *class;
    int nclasses;
    int nseeds;
    float *seed;
    double *weight;
    if ((o->seed_img_name == NULL) && (o->seed_txt_name == NULL)) 
      { /* No seed pixels are specified; mark all pixels as seeds: */
        if (o->verbose) { fprintf(stderr, "setting all seed classes to %d...\n", pnmclassif_DEFAULT_SEED_CLASS); }
        pnmclassif_set_all_classes(G, class, pnmclassif_DEFAULT_SEED_CLASS);
        /* Set root costs of seed pixels: */
        pnmclassif_set_root_weights_from_classes(G, rgb, input_chns, class, o->rfn);
      }
    else
      { /* Mark the seed pixels given as a PNM image: */
        if (o->seed_img_name != NULL)
          { pnmclassif_get_seeds_from_image_file(o, G, class); }
        else
          { pnmclassif_set_all_classes(G, class, pnmclassif_NON_SEED_CLASS); }
        pnmclassif_set_root_weights_from_classes(G, rgb, input_chns, class, o->rfn);
        /* Modify the seed pixels given as a text file: */
        if (o->seed_txt_name != NULL)
          { if (o->verbose) { fprintf(stderr, "reading seed file %s...\n", o->seed_txt_name); }
            pnmclassif_get_seeds_and_weights_from_text_file(o, G, rgb, input_chns, class, o->rfn);
          }
      }

    /* Compute the IFT: */
    if (o->verbose) { fprintf(stderr, "computing the IFT...\n"); }
    pnmclassif_compute_classification(G, rgb, input_chns, o->dfn, o->pfn, o->tbreak, o->order, &actual_max_cost);
      
    /* Write the output images: */
    
    /* Cost layer: */
    if (o->dist_img_name != NULL)
      { pnmclassif_write_cost_image(o, G, actual_max_cost); }
      
    /* Class layer: */
    if (o->class_img_name != NULL)
      { pnmclassif_write_class_image(o, G, class); }
    
    /* Predecessor layer: */
    if (o->prob_img_name != NULL)
      { pnmclassif_write_pred_image(o, G); }
      
    /* Roots indicator file: */
    if (o->root_img_name != NULL)
      { pnmclassif_write_root_image(o, G); }
      
    /* Root-spread image: */
    if (o->spread_img_name != NULL)
      { pnmclassif_write_spread_image(o, G, rgb, input_chns, input_maxval); }
      
    /* Single-class image: */
    if (o->extract_img_name != NULL)
      { pnmclassif_write_extract_image(o, G, class, rgb, input_chns, input_maxval); }

    /* Region bounding boxes: */
    if (o->class_moments_name != NULL)
      { pnmclassif_write_class_moments(o, G, class); }
      
    /* Postscript plot: */
    if (o->plot_eps_name != NULL)
      { pnmclassif_write_postscript_plot(o, G, rgb); }
      
    /* All done: */
    if (o->verbose) { fprintf(stderr, "done.\n"); }
    exit(0);
  }

void pnmclassif_compute_classification
  ( ift_graph_t *G, 
    frgb_t rgb[],
    int chns,
    pnmclassif_dist_fn_t *dfn, 
    pnmclassif_path_cost_fn_t *pfn, 
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
            double aC = dfn(rgb[si], a, rgb[ti], chns);
            return pfn(s->C, aC);
          }
      }
  }

void pnmclassif_read_pnm_input_image
  ( options_t *o, 
    uint16_t *input_maxvalP, 
    int *input_chnsP, 
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
    pnmclassif_set_values_from_image(input_img, G, rgb);
    (*input_maxvalP) = input_img->maxval;
    (*input_chnsP) = input_img->chns;
    (*GP) = G;
    (*rgbP) = rgb;
    uint16_image_free(input_img);
  }
    
void pnmclassif_set_all_classes(ift_graph_t *G, uint16_t class[], uint16_t lab)
  { uint64_t LL = lab; /* To supress a bogus "always false" gcc warning. */
    demand(LL <= pnmclassif_MAX_CLASS, "bad class");
    int32_t i;
    for (i = 0; i < G->nodes; i++) { class[i] = lab; }
  }
  
void pnmclassif_set_root_weights_from_classes
  ( ift_graph_t *G, 
    frgb_t rgb[], 
    int chns, 
    uint16_t class[],
    pnmclassif_root_cost_fn_t *rfn
  )
  { int32_t i;
    for (i = 0; i < G->nodes; i++) 
      { if (class[i] == pnmclassif_NON_SEED_CLASS)
          { G->node[i].C = +INF; }
        else
          { G->node[i].C = rfn(rgb[i], chns); }
      }
  }

void pnmclassif_get_seeds_from_image_file
  ( options_t *o, 
    ift_graph_t *G, 
    uint16_t class[]
  )
  { uint16_image_t *seed_img;
    seed_img = uint16_image_read_pnm_named(o->seed_img_name, o->verbose);
    if (seed_img->chns != 1) { pnm_error("seed image must be monochromatic"); }
    if (o->verbose) { fprintf(stderr, "marking seeds in graph...\n"); }
    pnmclassif_get_seeds_from_image(seed_img, G, class, o->verbose);
    uint16_image_free(seed_img);
  }

void pnmclassif_get_seeds_and_weights_from_text_file
  ( options_t *o, 
    ift_graph_t *G, 
    frgb_t rgb[], 
    int chns, 
    uint16_t class[],
    pnmclassif_root_cost_fn_t *rfn
  )
  {
    FILE *f = open_read(o->seed_txt_name, o->verbose);
    while (ch != EOF)
      { bool_t ok = fget_test_comment_or_eol(rd, '#', NULL);
        if (ok) { continue; }
        if (fget_test_eof(rd)) { break; }
        char ch = fget_char(rd);
        if (ch == '(')
          { int col = fget_int32(f); 
            demand((col >= 0) && (col < G->cols), "bad col in seed file"); 
            int row = fget_int32(f);
            demand((row >= 0) && (row < G->rows), "bad row in seed file");
            fget_skip_spaces_and_match(f, ")");
            fget_skip_spaces_and_match(f, "=");
            int lab = fget_int32(f);
            demand((lab >= 0) && (lab <= pnmclassif_MAX_CLASS), "bad class in seed file");
            int ig = ift_node_index(G, col, row);
            ift_node_t *pg = &(G->node[ig]);
            assert(pg->col == col);
            assert(pg->row == row);
            class[ig] = lab;
            /* Read handicap if any, else uses handicap = 0: */
            fget_skip_spaces(f);
            int ch1 = fgetc(f);
            if (ch1 != EOF) { ungetc(ch1, f); }
            if ((ch1 == '+') || (ch1 == '-') || (ch1 == '.') || ((ch1 >= '0') && (ch1 <= '9')))
              { pg->C = fget_double(f); }
            else
              { pg->C = rfn(rgb[ig], chns); }
            skip_rest_of_line(f);
            if (o->verbose)
              { fprintf(stderr, "  ( %5d %5d ) = %d %f\n", col, row, class[ig], pg->C); }
          }
        else
          { fprintf(stderr, "ch = «%c»\n", ch);
            demand(FALSE, "unexpected char in seed file");
          }
      }
    fclose(f);
  }

void pnmclassif_write_cost_image(options_t *o, ift_graph_t *G, ift_path_cost_t actual_max_cost)
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
    uint16_image_t *dist_img = pnmclassif_get_cost_image(G, o->max_cost_value, o->max_cost_pixel);
    uint16_image_write_pnm_named(o->dist_img_name, dist_img, 0, o->verbose);
    uint16_image_free(dist_img);
  }

void pnmclassif_write_class_image(options_t *o, ift_graph_t *G, uint16_t class[])
  { if (o->verbose) { fprintf(stderr, "extracting the class image...\n"); }
    /* Find maximum class (minimum 1): */
    uint16_t maxclass = 1;
    int i;
    for (i = 0; i < G->nodes; i++) { if (class[i] > maxclass) { maxclass = class[i]; } }
    uint16_image_t *class_img = pnmclassif_get_class_image(G, class, maxclass);
    uint16_image_write_pnm_named(o->class_img_name, class_img, 0, o->verbose);
    uint16_image_free(class_img);
  }

void pnmclassif_write_pred_image(options_t *o, ift_graph_t *G)
  { 
    if ((int)PNM_FILE_MAX_MAXVAL < 65535) 
      { fprintf(stderr, "PNM_FILE_MAX_MAXVAL = %d too small - pred image not generated\n", PNM_FILE_MAX_MAXVAL); }
    else
      { if (o->verbose) { fprintf(stderr, "extracting the pred image...\n"); }
        uint16_image_t *prob_img = pnmclassif_get_pred_image(G);
        uint16_image_write_pnm_named(o->prob_img_name, prob_img, 0, o->verbose);
        uint16_image_free(prob_img);
      }
  }

void pnmclassif_write_root_image(options_t *o, ift_graph_t *G)
  { if (o->verbose) { fprintf(stderr, "generating the root mask...\n"); }
    uint16_image_t *root_img = pnmclassif_get_root_image(G);
    uint16_image_write_pnm_named(o->root_img_name, root_img, 0, o->verbose);
    uint16_image_free(root_img);
  }

void pnmclassif_write_spread_image(options_t *o, ift_graph_t *G, frgb_t rgb[], int chns, uint16_t maxval)
  { uint16_image_t *spread_img;
    if (o->verbose) { fprintf(stderr, "extracting the spread image...\n"); }
    spread_img = pnmclassif_get_spread_image(G, rgb, chns, maxval);
    uint16_image_write_pnm_named(o->spread_img_name, spread_img, 0, o->verbose);
    uint16_image_free(spread_img);
  }

void pnmclassif_write_extract_image(options_t *o, ift_graph_t *G, uint16_t class[], frgb_t rgb[], int chns, uint16_t maxval)
  { if (o->verbose) { fprintf(stderr, "extracting the sub-image with class = %d...\n", o->extract_class); }
    uint16_image_t *extract_img = pnmclassif_get_single_class_image(G, class, o->extract_class, rgb, o->bgColor, chns, maxval);
    if (o->verbose) { fprintf(stderr, "maxval = %d\n", extract_img->maxval); }
    uint16_image_write_pnm_named(o->extract_img_name, extract_img, 0, o->verbose);
    uint16_image_free(extract_img);
  }

void pnmclassif_write_class_moments(options_t *o, ift_graph_t *G, uint16_t class[])
  { FILE *wr = open_write(o->class_moments_name, o->verbose);
    if (o->verbose) { fprintf(stderr, "writing the bounding boxes...\n"); }
    /* Find maximum seed class: */
    uint16_t maxclass = 0;
    int i;
    for (i = 0; i < G->nodes; i++) { if (class[i] > maxclass) { maxclass= class[i]; } }
    pnmclassif_write_boxes(wr, G, class, maxclass, o->box_margin);
    fclose(wr);
  }

void pnmclassif_write_postscript_plot(options_t *o, ift_graph_t *G, frgb_t rgb[])
  { double hsize = 12.0*(double)G->cols;
    double vsize = 12.0*(double)G->rows;
    if (o->verbose) { fprintf(stderr, "generating the tree plot %s...\n", o->plot_eps_name); }
    FILE *psfile = open_write(o->plot_eps_name, o->verbose);
    PSStream *ps = pswr_new_stream(NULL, psfile, TRUE, NULL, "letter", FALSE, hsize + 4.0, vsize + 4.0);
    pswr_new_canvas(ps, NULL);
    pswr_set_window
      ( ps, 
        0.0, (double)G->cols, 0.0, (double)G->rows,
        2.0, 2.0 + hsize, 2.0, 2.0 + vsize
      );
    pswr_set_grid(ps, G->cols, G->rows);
    ift_plot_pixel_values(ps, G, rgb, 0.5);
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.5, 0.0, 0.0); 
    ift_plot_forest_edges(ps, G);
    pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0.0, 0.0); 
    frgb_t root_rgb = (frgb_t){{ 1,1,1 }};
    ift_plot_forest_nodes(ps, G, TRUE,  1.0,  &root_rgb);
    frgb_t plain_rgb = (frgb_t){{ 0,0,0 }};
    ift_plot_forest_nodes(ps, G, FALSE, 0.5,  &plain_rgb);
    pswr_close_stream(ps);
  }

options_t *pnmclassif_parse_options(int argc, char **argv)
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
        o->rfn = pnmclassif_root_cost_fn_from_name(rfn_str);
      }
    else  
      { argparser_error(pp, "must specify a root cost function"); }
    
    if (argparser_keyword_present(pp, "-arcCost"))
      { char *dfn_str = argparser_get_next(pp);
        o->dfn = pnmclassif_dist_fn_from_name(dfn_str);
      }
    else  
      { argparser_error(pp, "must specify an arc cost function"); }
    
    if (argparser_keyword_present(pp, "-pathCost"))
      { char *pfn_str = argparser_get_next(pp);
        o->pfn = pnmclassif_path_cost_fn_from_name(pfn_str);
      }
    else  
      { argparser_error(pp, "must specify a path cost function"); }
    
    if (argparser_keyword_present(pp, "-output"))
      { /* Define matched filenames for costs, classes, and preds images: */
        char *output_str = argparser_get_next(pp);
        o->dist_img_name = txtcat(output_str, "-costs.pgm"); 
        o->class_img_name = txtcat(output_str, "-classes.pgm"); 
        o->prob_img_name = txtcat(output_str, "-preds.pgm"); 
      }
    else  
      { /* Allow individual filenames for costs, classes, and preds images: */
        if (argparser_keyword_present(pp, "-costs"))
          { o->dist_img_name = argparser_get_next(pp); }
        else  
          { o->dist_img_name = NULL; }

        if (argparser_keyword_present(pp, "-preds"))
          { o->prob_img_name = argparser_get_next(pp); }
        else  
          { o->prob_img_name = NULL; }

        if (argparser_keyword_present(pp, "-classes"))
          { o->class_img_name = argparser_get_next(pp); }
        else  
          { o->class_img_name = NULL; }
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
      { o->extract_class = argparser_get_next_int(pp, 0, INT_MAX);
        o->extract_img_name = argparser_get_next(pp);
      }
    else
      { o->extract_class = 0;
        o->extract_img_name = NULL; 
      }

    if (argparser_keyword_present(pp, "-boxes"))
      { o->class_moments_name = argparser_get_next(pp); }
    else  
      { o->class_moments_name = NULL; }
      
    if 
      ( (o->dist_img_name == NULL) &&
        (o->prob_img_name == NULL) &&
        (o->class_img_name == NULL) &&
        (o->root_img_name == NULL) &&
        (o->spread_img_name == NULL) &&
        (o->class_moments_name == NULL)
      )
      { fprintf(stderr, "what, no output files?\n"); }
    
    if (argparser_keyword_present(pp, "-margin"))
      { o->box_margin = argparser_get_next_int(pp, 0, INT_MAX); }
    else  
      { o->box_margin = 0; }
    
    if (argparser_keyword_present(pp, "-plot"))
      { o->plot_eps_name = argparser_get_next(pp); }
    else  
      { o->plot_eps_name = NULL; }

    if (argparser_keyword_present(pp, "-bgColor"))
      { int ch;
        for (ch = 0; ch < 3; ch++) { o->bgColor[ch] = argparser_get_next_int(pp, 0, INT_MAX); }
      }
    else if (argparser_keyword_present(pp, "-bgGray"))
      { uint16_t val = argparser_get_next_int(pp, 0, INT_MAX);
        int ch;
        for (ch = 0; ch < 3; ch++) { o->bgColor[ch] = val; }
      }
    else
      { int ch;
        for (ch = 0; ch < 3; ch++) { o->bgColor[ch] = 0; }
      }

    if (argparser_keyword_present(pp, "-maxCostValue"))
      { o->max_cost_value = argparser_get_next_int(pp, 0, INT_MAX); }
    else  
      { o->max_cost_value = 0; /* 0 means actual max */ }
    
    if (argparser_keyword_present(pp, "-maxCostPixel"))
      { o->max_cost_pixel = argparser_get_next_int(pp, 1, PNM_FILE_MAX_MAXVAL); }
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
      { int r; char c;
        do { r = fgetc(f); c = (char)r } while ((r != EOF) || (c != '\n')); }
    else
      { fget_eol(f); }
  }
