/* Test of jspng.h, uint16_image_io_png.h */
/* Last edited on 2017-06-25 16:25:34 by stolfilocal */

#define PROG_NAME "pnmtoobj"
#define PROG_DESC "converts a height map and/or a color image to a Wavefront OBJ file"
#define PROG_VERS "1.0"

#define pnmtoobj_C_COPYRIGHT \
  "Copyright © " pnmtoobj " by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -geometry {G_IMAGE} ] \\\n" \
  "    [ -texture {T_IMAGE} [ -gammma {GAMMA} {BIAS} ] ] \\\n" \
  "    [ -scale {SCX} {SCY} {SCZ} ] \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  The program reads a geometry map image {G_IMAGE} (such as a digital" \
  " terrain or a compacted quadrilateral mesh) and/or a " \
  " color (\"texture\") image {T_IMAGE}, all in a Netpbm" \
  " file format (\".pbm\", \".pgm\", or \".ppm\").  It writes them to standard" \
  " output as a 3D triangular mesh in the Wavefront \".obj\" format.\n" 
  "\n" \
  "  Basically, each pixel of the input image(s) is converted to four triangles of" \
  " the mesh.  The coordinates of the vertices of the mesh are taken" \
  " from {G_IMAGE}, or from a regular grid if that file is not" \
  " given.  The colors of the triangles are taken from" \
  " the {T_IMAGE}, or set to a uniform grayish tone if that file not given.\n" \
  "\n" \
  "  The program ignores the extensions of the input filenames, and" \
  " interprets the contents according to the 16-bit \"magic number\" at" \
  " the start of the file.  All sample values are converted from integers" \
  " in {0..MAXVAL} to floating-point values" \
  " between 0 and 1; where {MAXVAL} is the maximum sample value as" \
  " declared in the input file.   Note that in PBM (binary) files" \
  " the encoding is reversed (0-->1, 1-->0).\n" \
  "\n" \
  "OUTPUT MESH SIZE AND TOPOLOGY\n" \
  "\n" \
  "  The output mesh is a grid of /quadrangles/" \
  " with {NX} rows by {NY} columns.  Each quadrangle consists of four" \
  " triangles that share a common /central vertex/.\n" \
  "\n" \
  "  If the {G_IMAGE} is given, then each corner vertex corresponds to" \
  " a pixel of {G_IMAGE}.  Thus, if {G_IMAGE} has {GNX} columns" \
  " and {GNY} rows, the output mesh will have {NX=GNX-1} and {NY=GNY-1}.\n" \
  "\n" \
  "  If both the {G_IMAGE} and the {T_IMAGE} are given, the sizes" \
  " must be the same, or the texture must be smaller by one pixel" \
  " than the {G_IMAGE}, in both directions.  That is, if {T_IMAGE} has" \
  " {TNX} columns and {TNY} rows, then either {TNX=GNX} and {TNY=GNY}, or" \
  " {TNX=GNX-1} and {TNY=GNY-1}.  In either case, the sizes {NX,NY} of" \
  " the mesh will be {GNX-1,GNY-1}, as above.\n" \
  "\n" \
  "  If the {G_IMAGE} is not given, then each quadrangle of the mesh" \
  " corresponds to one pixel of the {T_IMAGE}, and the output mesh" \
  " will have {NX=TNX} and {NY=TNY}.\n" \
  "\n" \
  "OUTPUT MESH COORDINATES\n" \
  "\n" \
  "  The {G_IMAGE} is a 3-channel color image (PPM format, the red, green," \
  " and blue channel samples become the {X}, {Y}, and {Z} coordinates of" \
  " the corner vertex.\n" \
  "\n" \
  "  If the {G_IMAGE} is a single-channel grayscale image (PGM or PBM format)," \
  " the pixel samples define the {Z} coordinate of each corner vertex," \
  " while the {X} and {Y} coordinates are set to successive integers," \
  " respectively from 0 to {NX} and 0 to {NY}, inclusive.\n" \
  "\n" \
  "  In both cases, the integer sample values in the {G_IMAGE} are converted to" \
  " floating-point values by simple linear interpolation" \
  " from {0..MAXVAL} to {{0 _ 1]}, with no gamma correction.\n" \
  "\n" \
  "  If the {G_IMAGE} is not given, the {X} and {Y} coordinates of the" \
  " corner vertices are set to a regular grid, as above, and the {Z}" \
  " coordinates are all set to zero.\n" \
  "\n" \
  "  In all cases, the center vertex of each quadrangle will be the" \
  " average of the four corner vertices, whether these came from" \
  " the {G_IMAGE} or from the regular grid.\n" \
  "\n" \
  "OUTPUT MESH COLORS\n" \
  "  If the {T_IMAGE} is not given, the output mesh will have no" \
  " color information.\n" \
  "\n" \
  "  If the {T_IMAGE} is given, but the {G_IMAGE} is not given or" \
  " is 1 pixel larger than the {T_IMAGE} in {X} and {Y} , then each" \
  " pixel of the {T_IMAGE} defines the color of the four triangles" \
  " in the corresponding quadrangle of the mesh.  The colors then" \
  " will be discontinuous at the boundaries of the quadrangles.\n" \
  "\n" \
  "  If both the {T_IMAGE} and the {G_IMAGE} are given, and have the" \
  " same {X} and {Y} sizes, then each pixel of the {T_IMAGE} defines" \
  " the color of the corresponding corner vertex of the mesh.  The" \
  " color of each center vertex will be the average of the colors" \
  " of the four corner vertices.  Each triangle then should be" \
  " rendered by affine interpolation of its corner colors (Gouraud shading).\n" \
  "\n" \
  "  In any case, the integer sample values in the {T_IMAGE} are converted\n" \
  " to floating-point values by affine interpolation" \
  " from {0..MAXVAL} to {[EPS _ 1-EPS]}, where {EPS=0.5/(MAXVAL+1)}," \
  " followed by gamma correction.\n" \
  "\n" \
  "OPTIONS\n" \
  "  -geometry {G_IMAGE}\n" \
  "    This optional parameter is the name of the image file that defines" \
  " the geometry of the output mesh.  At least one" \
  " of \"-geometry\" and \"-texture\" must be specified.\n" \
  "\n" \
  "  -texture {T_IMAGE}\n" \
  "    This optional parameter is the name of the image file that" \
  " defines the colors of the mesh.  At least one" \
  " of \"-geometry\" and \"-texture\" must be specified.\n" \
  "\n" \
  "  -gamma {GAMMA} {BIAS}\n" \
  "    This optional parameter specifies the parameters of the" \
  " non-linear (\"gamma\") encoding of luminosity values in the texture image.   If omitted, the program assumes \"-gamma 2.2222 0.0327, which is a good approximation to the ITU-R BT.709 standard used by \n" \
  "\n" \
  "  -scale {SCX} {SCY} {SCZ}\n" \
  "    This optional parameter specifies a scale factor that can" \
  " be applied to each coordinate, after it has been defined" \
  " as above.  For a digital terrain  (grayscale {G_IMAGE}), for" \
  " example, \"-scale 1 1 10\" would set the maximum height of the terra  \n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "   pgmtopnm(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2017-03-03 by Jorge Stolfi, IC-UNICAMP, from bits of {fni_view}.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2017-03-03 J.Stolfi: created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " pnmtoobj_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <r3.h>
#include <float_image.h>
#include <float_image_io_pnm.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { char *geomety; /* Name of geometry image file. */
    char *texture; /* Name of texture (color) image file. */
    r3_t scale;    /* Coordinate scaling factors. */
  } options_t;

/* INTERNAL PROTOTYPES */

options_t *pnmtoobj_parse_options(int argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int main(int argc,char** argv);

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  {
    /* Parse command line: */
    options_t *o = pnmtoobj_parse_options(argc, argv);
    
    /* Read geometry and texture images: */
    bool_t isMask = FALSE;
    float_image_t *gimg = NULL; /* Geometry image. */
    float_image_t *timg = NULL; /* Texture image. */
    if (o->geometry != NULL) 
      { double gamma = 1.0; /* Gamma correction exponent. */
        double bias = 0.0;    /* Bias for gamma correction. */
        gimg = float_image_read_pnm_named(o->geometry, isMask, gamma, bias, TRUE, TRUE, FALSE); 
      }
    if (o->etxture != NULL) 
      { double gamma = 1/0.4500; /* Gamma correction exponent. */
        double bias = 0.0327;    /* Bias for gamma correction. */
        timg = float_image_read_pnm_named(o->texture, isMask, gamma, bias, TRUE, TRUE, FALSE);
      }
    
    float_image_t *timg = float_image_read_pnm_named(o->, isMask, gamma, bias, TRUE, TRUE, FALSE);
    
    r3_t *vert = NULL; /* The vertices of the mesh. */
    i3_t *face;        /* The triangles. */
    frgb_t *fcol,      /* The colors of triangles. */
    void pto_convert_height_map(ht, c, zscale, tx, &vert, &face, &tint);

    /* Modify it: */
    fprintf(stderr, "modifying image...\n");
    uint16_image_t *omg = uint16_image_new(img->cols, img->rows, img->chns);
    frobnicate_image(img, omg);
    fprintf(stderr, "------------------------------------------------------------\n");
    uint16_image_describe(stderr, "modified", omg);
    fprintf(stderr, "------------------------------------------------------------\n");

    /* Write it out: */
    double oGamma = iGamma;
    uint16_image_io_png_write(oName, omg, oGamma, TRUE);
    
    /* Read it back: */
    double bGamma;
    uint16_image_t *bmg = uint16_image_io_png_read(oName, &bGamma, TRUE);
    
    /* Check attributes: */
    fprintf(stderr, "------------------------------------------------------------\n");
    uint16_image_describe(stderr, oName, bmg);
    fprintf(stderr, "------------------------------------------------------------\n");
    compare_images(omg, bmg);
    if (isnan(oGamma))
      { demand(isnan(bGamma), "gamma should be {NAN}"); }
    else
      { demand(fabs(bGamma-oGamma)/sqrt(bGamma*oGamma) < 0.000001, "inconsistent {gamma}"); }
    
    uint16_image_free(img);
    uint16_image_free(omg);
    uint16_image_free(bmg);
    
    free(iName);
    free(oName);
    fprintf(stderr, "test successful!\n");
    fprintf(stderr, "============================================================\n");
    
    return 0;
  }

void pto_convert_height_map
  ( float_image_t *ht, 
    int c, 
    double zscale,
    float_image_t *tx,
    r3_t **vert,
    i3_t **face,
    frgb_t **tint
  )
  {
    /* Get grid sizes: */
    demand((ht != NULL) || (tx != NULL), "either height map or texture map must be given");
    int HNC, HNX, HNY;  /* Height map sizes. */
    int TNC, TNX, TNY;  /* Texture map sizes. */
    int NC, NX, NY; /* Triangle grid sizes. */
    if (ht != NULL) 
      { float_image_get_size(ht, &HNC, &HNX, &HNY);
        demand(HNC == 1, "Height map must be grayscale");
        /* Triangle grid size is that of the terrain: */
        NC = HNC; NX = HNX; NY = HNY;
      }
    if (tx != NULL) 
      { float_image_get_size(tx, &TNC, &TNX, &TNY);
        if (hx == NULL) 
          { /* Triangle grid size is that of texture image: */
            NC = TNC; NX = TNX; NY = TNY;
          }
        else
          { demand((TNX == NX-1) && (TNY == NY-1), "texture and heigh map size mismatch"); }
      }
    else
      { /* Pretend texture has the right size: */
        TNC = 3; TNX = NX-1; TNY = NY-1;
      }

    /* We need buffers for two rows of height map samples: */
    float va[NX], vb[NX];
    float *v0 = va; /* Current row of height samples. */
    float *v1 = vb; /* Next row of height samples. */
    
    /* We need a buffer for one row of texture map colors: */
    float clr[TNC*TNX];

    /* Get first row of height samples (or all zero if no height map): */
    if (ht != NULL)
      { float_image_get_sample_row(ht, c, 0, HNX-1, 0, v0); }
    else
      { for(int x; x < NX; x++) { v0[x] = 0.0; } }
    
    /* Scan rows of triangle array: */
    int x, y;
    for(y = 0; y < NY-1; y++)
      { /* Get next row of height samples (or all zero if no height map): */
        if (ht != NULL)
          { float_image_get_sample_row(ht, c, 0, HNX-1, y+1, v1);}
        else
          { for(int x; x < NX; x++) { v1[x] = 0.0; } }
        
        /* Get next row of colors, or all greyish if no image: */
        if (tx != NULL)
          { float_image_get_pixel_row(tx, 0, TNX-1, y, clr); }
        else
          { for (int x = 0; x < TNX; x++) 
              { for (int c = 0; c < TNC; c++) 
                  { clr[x*TNC + 0] = 0.900f; 
                    clr[x*TNC + 1] = 1.000f; 
                    clr[x*TNC + 2] = 0.950f; 
                  }
              }
          }
        
        /* Now convrt the pixels: */
        float *pv0 = v0; /* Pointer to height in col {x}, row {y}. */
        float *pv1 = v1; /* Pointer to height in col {x}, row {y+1}. */
        float *pc = clr; /* Pointer to color of pixel in col {x}, row {y}. */
        for(x = 0; x < NX-1; x++)
          { /* Get color in GL format: */
            GLfloat CR, CG, CB;
            if (TNC == 1)
              { CR = CG = CB = pc[0]; }
            else
              { CR = pc[0]; CG = pc[1]; CB = pc[2]; }
            pc += TNC;
            /* Get heights at four pixel corners in GL format: */
            float z00 = (float)(zscale*pv0[0]);  /* Height at {(x,   y  )}. */
            float z10 = (float)(zscale*pv0[1]);  /* Height at {(x+1, y  )}. */
            float z01 = (float)(zscale*pv1[0]);  /* Height at {(x,   y+1)}. */
            float z11 = (float)(zscale*pv1[1]);  /* Height at {(x+1, y+1)}. */
            pto_convert_cell(x,y, z00,z10,z01,z11, CR,CG,CB);
            pv0++; pv1++;
          }
        /* Swap row buffers, {v1} now becomes {v0}: */
        { float *t = v0; v0 = v1; v1 = t; }
      }
  }
            
void pto_convert_cell
  ( FILE *wr,
    int x, int y, 
    float z00, float z10, float z01, float z11, 
    float CR, float CG, float CB
  )
  {
    /* Low and high coordinates of pixel on XY plane: */
    float x0 = (float)x;
    float y0 = (float)y;
    float x1 = x0 + 1.00f;
    float y1 = y0 + 1.00f;
    
    fprintf(wr, "v %8.2f %8.2f %8.2f\n", 
    /*
      Plot pixel as four triangles {p,q,u} where {u} is the center 
      and {p,q} scan the boundary in the order
        01<----11
        |       ^
        |  .u   |
        V       |
        00---->10
    */

    /* The height at the center is the average of the four corner heights: */
    float xu = (float)(x + 0.5);
    float yu = (float)(y + 0.5);
    float zu = (z00 + z10 + z01 + z11)/4;

    /* !!! Should use {GL_TRIANGLE_FAN} instead of {GL_TRIANGLES}. !!! */
    glBegin(GL_TRIANGLES);
    pto_output_triangle(xu, yu, zu, x0, y0, z00, x1, y0, z10);
    pto_output_triangle(xu, yu, zu, x1, y0, z10, x1, y1, z11);
    pto_output_triangle(xu, yu, zu, x1, y1, z11, x0, y1, z01);
    pto_output_triangle(xu, yu, zu, x0, y1, z01, x0, y0, z00);
    glEnd();
  }

options_t *pnmtoobj_parse_options(int argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the command line argument record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    /* Parse keyword parameters: */
    
    /* Parse the boolean option {op1}: */
    o->op1 = argparser_keyword_present(pp, "-op1");
    
    /* Parse the string option {op2}: */
    if (argparser_keyword_present(pp, "-op2"))
      { o->op2 = argparser_get_next(pp); }
    else
      { o->op2 = "NONE"; }

    /* Parse positional arguments: */
    argparser_skip_parsed(pp);
    
    if (argparser_next(pp) != NULL)
      { o->infile = argparser_get_next(pp); }
    else
      { o->infile = "-"; }

    if (argparser_next(pp) != NULL)
      { o->outfile = argparser_get_next(pp); }
    else
      { o->outfile = "-"; }

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

