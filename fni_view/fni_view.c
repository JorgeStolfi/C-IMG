#define PROG_NAME "fni_view"
#define PROG_DESC "3D interactive visualization of a height map"
#define PROG_VERS "2.0"

#define PROG_C_COPYRIGHT "Copyright © 2005 Universidade Estadual Fluminense (UFF)."

/* Last edited on 2024-12-23 09:09:26 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -range {VMIN} {VMAX} | -range auto ] \\\n" \
  "    [ -scale {VSCALE} ] \\\n" \
  "    [ -isMask {ISMASK} ] \\\n" \
  "    [ -colorize {VMIN} {VMAX} | -colorize auto ] \\\n" \
  "    [ -channel {CHAN} ] \\\n" \
  "    [ -texture {TEXTURE_FILE} | -hist {HISTFLAG}  ] \\\n" \
  "    [ < ] {HEIGHT_FILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Reads an image from file {HEIGHT_FILE} and interprets it as a" \
  " height map.  Displays a selected channel of it" \
  " as a terrain in 3D or as a 3D bar histogram.\n" \
  "\n" \
  "  If \"-texture\" is specified, reads also a second image from" \
  " {TEXTURE_FILE}, and uses it asa texture map to colorize the terrain.\n" \
  "\n" \
  "  The interface allows the user to move the camera around with the" \
  " mouse, and modify the view by typing single-charater commands to the window.\n" \
  "\n" \
  "  The sample value in channel {CHAN}, column {x}, row {y} of the input image is" \
  " taken as the Z coordinate of a 3D terrain, above the" \
  " point {(x,y)} of the {Z=0} plane.  Normally the terrain is" \
  " displayed by conecting those points with quadrangular" \
  " patches formed by 4 or 8 triangles.  Alternatively, they may be" \
  " displayed as parallelepipeds (bar-graph style). The sample" \
  " values may be optionally magnified" \
  " by a scale factor before displaying.\n" \
  "\n" \
  "  The texture image may have either the same dimensions as the" \
  " height image, or exactly one column and one row less than it.  In the" \
  " first case, the texture pixel in column {x}, row {y} is taken" \
  " as the color of the terrain above the unit square with" \
  " diagonal {(x-1/2,y-1/2)} and {(x+1/2,y+1/2)}.  In the" \
  " second case, the pixel in column {x}, row {y} of the texture map is" \
  " taken as the color of the unit-side square whose diagonal" \
  " corners are {(x,y)} and {(x+1,y+1)}.\n" \
  "\n" \
  " If no texture map is specified, the terrain is colorized" \
  " with colors computed from the sample values.\n" \
  "\n" \
  "  The input files may be in the PNM (PBM/PGM/PPM) format" \
  " (see {uint16_image.h}) or in the  float-valued multi-channel" \
  " image format (FNI) (see {float_image.h}).  If the" \
  " argument {HEIGHT_FILE} is \"-\" or omitted, the program reads it" \
  " from standard input, assuming the FNI format.\n" \
  "\n" \
  "  If the main image is in the PNM format, its samples are" \
  " converted to floats in the range {[0_1]} as specified by" \
  " the \"-isMask\" option. \n" \
  "\n" \
  "OPTIONS\n" \
  "  -channel {CHAN}\n" \
  "    Specifies the channel of the height file that is to be used as the" \
  " height map.  Channel numbers start at 0; the program fails if the" \
  " channel does not exist in the input file.  The default is \"-channel 0\".  See" \
  " also the interactive keyboard commands 'c' and '0' to '9'.\n" \
  "\n" \
  "  -isMask {ISMASK}\n" \
  "    This optional Boolean argument specifies the interpretation" \
  " of integer sample values when the main input file is in the PNM format, specifically" \
  " how the integer values {0..MAXVAL} are mapped to the range {[0_1]}.  If {ISMASK} is true (\"T\" or 1)," \
  " " sample_conv_0_1_isMask_true_INFO "  If {ISMASK} is false (\"F\" or 0)," \
  " " sample_conv_0_1_isMask_false_INFO "  The default is \"-isMask F\". This" \
  " option does not affect the reading of PNM texture files.\n" \
  "\n" \
  "  -texture {TEXTURE_FILE}\n" \
  "    Specifies the name of the FNI file containing the texture" \
  " map.  If {TEXTURE_FILE} is \"-\", reads the texture map from" \
  " standard input too (after the height map, if applicable).  See" \
  " also the interactive keyboard command 't'.\n" \
  "\n" \
  "    The texture map should have one column and one row" \
  " less than the height map.\n" \
  "\n" \
  "  -hist {HISTFLAG}\n" \
  "    If {HISTFLAG} is true (\"T\" or 1), the map is to be drawn" \
  " as an histogram, with each pixel as a horizontal square and vertical" \
  " walls between adjacent pixels. In this case the \"-texture\" option" \
  " is ignored. If {HISTFLAG} is false (\"F\" or 0), draws values as a continuous" \
  " triangle mesh.  The default is \"-hist F\".  See" \
  " also the interactive keyboard command 'h'.\n" \
  "\n" \
  "  -range {VMIN} {VMAX}\n" \
  "  -range auto\n" \
  "    Specifies the nominal min and max pixel values for the purpose of" \
  " scaling and centering.  In the first variant, the heights" \
  " obtained from the height map will be scaled so that the given range" \
  " from {XMIN} to {XMAX} will be stretched to be equal to the diameter" \
  " of the map's domain.  The \"auto\" option is similar but" \
  " assumes {VMIN} and {VMAX} are the minimum and maximum" \
  " height values.  The default is \"-range auto\".  See also" \
  " the interactive keyboard commands 's' and 'S'.  This argument" \
  " also sets the position of the" \
  " lower and upper reference planes; see the 'p' key command.\n" \
  "\n" \
  "  -scale {VSCALE}\n" \
  "    Specifies an additional multiplying scale factor for the heights" \
  " recovered from the height map.  This scaling factor is applied in" \
  " addition to the scaling implied \"-range\", if given.  The default" \
  " is \"-scale 1.0\".  See" \
  " also the interactive keyboard commands 's' and 'S'.\n" \
  "\n" \
  "  -colorize {VMIN} {VMAX}\n" \
  "  -colorize auto\n" \
  "    Specifies the range of sample values to assume when computing" \
  " pseudocolors.  This argument is used only if the texture map is" \
  " omitted or turned off.  In the first option, samples" \
  " in the range {[VMIN _ VMAX]} are mapped to colors of an appropriate" \
  " pseudocolor scale; samples outside that range are mapped to the end" \
  " points of the scale.  The \"auto\" option will set {VMIN} and {VMAX} to" \
  " the maximum and minimum" \
  " samples in the current channel.  The default is \"-colorize auto\".  The" \
  " color mapping is not affected by the \"-range\" and \"-scale\" arguments" \
  " or the 's' and 'S' command keys.\n" \
  "\n" \
  "VIEWING INTERACTION\n" \
  "  The viewing direction can be changed by dragging with the" \
  " left mouse button.  The following keyboard commands" \
  " also are recognized:\n" \
  "\n" \
  "  'q' : Quit.\n" \
  "  'Q' : Quit.\n" \
  "  ESC : Quit.\n" \
  "  't' : Toggle texturization.\n" \
  "  'p' : Toggle display of the reference plane.\n" \
  "  'o' : Toggle display of zero-height plane.\n" \
  "  'h' : Toggle between histogram and terrain modes.\n" \
  "  'Z' : Zoom in, down to a limit.\n" \
  "  'z' : Zoom out, up to a limit.\n" \
  "  'S' : Increase the height scale factor.\n" \
  "  's' : Reduce height scale factor.\n" \
  "  'c' : Cycle through channels of the image.\n" \
  "  '0' to '9': Display the indicated channel of the image.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnm_to_fni(1), fni_to_pnm(1), {float_image.h}.\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2005 by Rafael Saracchini, IC-UFF.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  apr/2006 by Jorge Stolfi, IC-UNICAMP: Recast to use jslibs, argparser, etc..\n" \
  "  jul/2009 by R. Saracchini: Changed to use lights and display list.\n" \
  "  jun/2010 by R. Saracchini: Added mouse support.\n" \
  "  jul/2010 by Jorge Stolfi: Added PGM/PPM input.\n" \
  "  jul/2010 by Jorge Stolfi: Added support for same-size textures.\n" \
  "  jul/2010 by Jorge Stolfi: Added \"-colorize\" option.\n" \
  "  aug/2010 by Jorge Stolfi: Added \"-isMask\" option.\n" \
  "  jan/2011 by Jorge Stolfi: Added \"-hist\" option.\n" \
  "  oct/2022 by Jorge Stolfi: Added \"-range\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <assert.h>
#include <string.h>

#define __APPLE__ 0
#include <GL/glu.h>
#include <GL/glut.h>

#include <argparser.h>
#include <float_image.h>
#include <float_image_read_pnm.h>
#include <sample_conv.h>
#include <frgb.h>
#include <frgb_path.h>
#include <bool.h>
#include <jsfile.h>
#include <r3.h>

#include <fvw_GL.h>
#include <fvw_paint.h>
#include <fvw_paint_node_colored.h>
#include <fvw_paint_cell_colored.h>
#include <fvw_paint_self_colored.h>
#include <fvw_paint_self_colored_hist.h>

#define INF INFINITY

#define fvw_zoom_step_ratio (15.0/16.0)
  /* Zoom in/out kbd commands multiply/divide the obs distance by this amount. */

#define fvw_ht_scale_step_ratio (pow(2.0,1.0/12.0))
  /* Scale incr/decr kbd commands multiply/divide the extra height scale factor by this amount. */

/* COMMAND-LINE OPTIONS */

typedef struct options_range_t 
  { bool_t vauto;   /* TRUE to get height range from height map values. */
    float vmin;     /* Min user-specified sample value, or {NAN}. */
    float vmax;     /* Max user-specified sample value, or {NAN}. */
  } options_range_t;
  /* User options for heights caling or colormapping pixel ranges. */
 
options_range_t fvw_parse_options_range(argparser_t *pp, char *keyword);

typedef struct options_t
  { char* height_file;        /* Name of height map. */
    /* Height scaling options: */
    options_range_t range;    /* Height range options for scaling. */
    double scale;             /* Additional height scale factor. */
    /* Colorizing options: */
    char* texture_file;       /* Name of texture map, or NULL if not specified. */
    options_range_t colorize; /* Height range options for color mapping. */
    /* Other options: */
    uint32_t channel;          /* Initial channel to display. */
    bool_t isMask;        /* TRUE if input samples are to be converted as in a mask. */
    bool_t hist;          /* TRUE specifies histogram-like plot. */
  } options_t;

/* GLOBAL VARIABLES USED BY GL METHODS */

typedef struct fvw_state_range_t
  { float def_vmin;        /* Default {vmin} for new channels or {NAN} if auto. */
    float def_vmax;        /* Default {vmax} for new channels or {NAN} if auto. */
    float_vec_t vmin;      /* High pix values per channel or {NAN} if undef (mutable). */
    float_vec_t vmax;      /* Low pix values per channel or {NAN} if undef (mutable). */
  } fvw_state_range_t;
  /* Specifies range parameters for height scaling or colormapping. */

fvw_state_range_t fvw_state_range_make(options_range_t *orn, int32_t NC);
  /* Creates a {fvw_state_range_t} from the user options {orn} for
    a height map image with {NC} channels. */
  
void fvw_get_vrange(uint32_t c, fvw_state_range_t *srn, float_image_t *ht, float *vminP, float *vmaxP);  
  /* Obtains the height scaling or colormapping pixel range {*vminP} and {*vmaxP} for channel {c}
    of the image {ht}.  If the range has been previously computed and saved in the record {srn},
    takes it from there.  Otherwise, if {srn} has a user-specified range, uses that range.
    Otherwise recomputes the range from the pixel values of {ht}. Either way, saves these
    values in the {srn} record for future use. */

typedef struct fvw_state_t 
  { 
    /* Height map: */
    float_image_t *ht;  /* Height image (const). */
    double RAD;         /* Nominal radius of terrain (for default height scaling). */

    /* Channel to display: */
    uint32_t channel;         /* Current channel being displayed (mutable). */

    /* Display style: */
    bool_t hist;         /* TRUE specifies histogram-like plot (mutable). */

    /* Reference pixel values for height scaling: */
    fvw_state_range_t ht_range; 

    /* Extra height magnification factor: */
    double def_ht_mag;        /* Default extra height scale factor for new channels (mutable). */
    double_vec_t ht_mag;      /* Current extra scale factor per channel, or NAN if undefined (mutable). */
    
    /* Texture map: */
    float_image_t *tx;   /* Texture image (const). */
    bool_t node_tx;      /* TRUE if {tx} is suitable for node-painting (const). */
    bool_t cell_tx;      /* TRUE if {tx} is suitable for cell-painting (const). */
    bool_t texturize;    /* When TRUE, use texture map (mutable). */
    
    /* Reference pixel values for colormapping: */
    fvw_state_range_t cm_range; 
   
    /* Reference plane: */
    bool_t refPlane;     /* When TRUE, show {ref_vmin} and {ref_vmax} planes (mutable). */
    bool_t orgPlane;     /* When TRUE, show reference plane placed at z = 0 (mutable). */
    
    /* Current window dimensions (set by {fvw_reshape_method}): */
    GLint window_HSize; /* Width (mutable). */
    GLint window_VSize; /* Height (mutable). */

    /* Observer's position relative to center of bbox: */
    GLfloat azimuth;    /* Azimuth from X axis (degrees, mutable). */
    GLfloat elevation;  /* Elevation from XY plane (degrees, mutable). */
    GLfloat distance;   /* Distance (pixels, mutable). */
    
    /* Mouse state: */
    int32_t mouse_x,mouse_y; /* Position of last mouse event (mutable). */
  } fvw_state_t;

static fvw_state_t *fvw_state = NULL; /* The state displayed in the GL window. */

/* INTERNAL PROTOTYPES */

fvw_state_t *fvw_create_state(options_t *o);
  /* Creates and initializes a window data record, suitable for use
    used by the methods {fvw_display_method} etc. */

float_image_t *fvw_read_image(char *name, bool_t isMask, int32_t *NCP, int32_t *NXP, int32_t *NYP);
  /* Reads a float image, in the format of {float_image_write} (if the
    {name} ends in ".fni") or in PBM/PGM/PPM format (if {name} ends in
    ".pbm", ".pgm", ".ppm", ".pnm"). If {name} is "-", reads from
    {stdin} assuming {float_image_write} format. The procedure also
    sets {*NCP,*NXP,*NYP} to the number of channels, columns, and rows.*/

float_image_t *fvw_read_float_image(char *name);
  /* Reads a float image, in the format of {float_image_write}, from
    the named file (or {stdin} if {name} is "-"). */

options_t *fvw_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int32_t main(int32_t argc,char** argv);

/* MAJOR PAINTING FUNCTIONS */

void fvw_paint_everything(fvw_state_t *w);
  /* Repaints the terrain according to the data and parameters in {w}. */

void fvw_set_perspective(fvw_state_t *w);
  /* Sets the GL perspective view parameters according to the 
    window size, observer position, and bounding box data 
    stored in {w}. */

void fvw_set_lights(fvw_state_t *w);
  /* Sets the GL lighting parameters according to the state {w}. */

/* MEDIUM-LEVEL PAINTING */

void fvw_paint_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double ht_scale, 
    float cm_vmin,
    float cm_vmax,
    bool_t hist,
    float_image_t *tx,
    bool_t node_tx,
    bool_t cell_tx
  );
  /* Paints the height map consisting of channel {c} of image {ht},
    with the Z coordinate being the sample values scaled by {ht_scale}.
    If {tx} is not null, it must have either the same size as {ht}, or
    one pixel less in both axes; in either case, paints the terrain
    with color {tx}. If {tx} is NULL, paints with colors based on the
    unscaled sample values relative to the unscaled sample range {[cm_vmin _ cm_vmax]}.
    If {hist} is true, paints as a histogram (and ignores {tx}). */

void fvw_paint_reference_plane(fvw_state_t *w, float z);
  /* Paints a horizontal rectangle at height {z} 
    (assumed to be already scaled). */

/* GL WINDOW METHODS */

void fvw_display_method(void);
  /* Processes redisplay requests by the window manager. */

void fvw_reshape_method(int32_t width, int32_t height);
  /* Processes window reshape events. */

void fvw_keyboard_method(unsigned char key, int32_t x, int32_t y);
  /* Processes keyboard events. */
     
void fvw_passivemouse_method( int32_t x, int32_t y);
  /* Processes passive mouse events (no mouse clicked) */

void fvw_activemouse_method( int32_t x, int32_t y);
  /* Processes active mouse events (mouse clicked) */

void fvw_special_method(int32_t key, int32_t x, int32_t y);
  /* Processes special user input events. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char** argv)
  {
    fvw_GL_initialize_libraries(&argc, argv);
  
    options_t *o = fvw_parse_options(argc, argv);

    fvw_state = fvw_create_state(o);
    
    fvw_GL_initialize_window
      ( o->height_file,
        fvw_display_method,
        fvw_reshape_method,
        fvw_keyboard_method,
	fvw_passivemouse_method,
	fvw_activemouse_method,
        fvw_special_method
      );
    
    fvw_GL_start_viewing();
    
    return 0;
  }

fvw_state_t *fvw_create_state(options_t *o)
  { 
    fvw_state_t *w = (fvw_state_t *)notnull(malloc(sizeof(fvw_state_t)), "no mem");

    w->window_HSize = fvw_GL_default_window_HSize;
    w->window_VSize = fvw_GL_default_window_VSize;

    /* Read height map, check and save the channel index: */
    int32_t ht_NC, ht_NX, ht_NY;
    w->ht = fvw_read_image(o->height_file, o->isMask, &ht_NC, &ht_NX, &ht_NY);
    fprintf(stderr, "height map has %d channels, %d columns, %d rows\n", ht_NC, ht_NX, ht_NY); 
    demand(ht_NC > 0, "invalid num of channels in height map");
    
    /* The initial channel is user-specified: */
    w->channel = o->channel;
    if ((w->channel < 0) || (w->channel >= ht_NC))
      { fprintf(stderr, "invalid height map channel %d, using 0\n", w->channel);
        w->channel = 0;
      }

    /* Pixel ranges for height scaling: */
    w->ht_range = fvw_state_range_make(&(o->range), ht_NC);

    /* Extra height scaling factors: */
    w->def_ht_mag = o->scale;
    w->ht_mag = double_vec_new((uint32_t)ht_NC);
    for (uint32_t c = 0;  c < ht_NC; c++) { w->ht_mag.e[c] = NAN; }

    /* Pixel ranges for color mapping: */
    w->cm_range = fvw_state_range_make(&(o->colorize), ht_NC);
    
    /* The histogram-style flag is initially obtained from the command line: */
    w->hist = o->hist;
     
    /* Read texture map, if any,and check its size: */
    int32_t tx_NC, tx_NX, tx_NY; /* Channels, columns and rows of texture map. */
    if (o->texture_file != NULL)
      { bool_t tx_isMask = FALSE; /* Assume that texture has smooth pixel distr. */
        w->tx = fvw_read_image(o->texture_file, tx_isMask, &tx_NC, &tx_NX, &tx_NY);
        fprintf(stderr, "texture map has %d channels, %d columns, %d rows\n", tx_NC, tx_NX, tx_NY); 
        demand((tx_NC == 1) || (tx_NC == 3), "invalid num of channels in texture map");
        w->node_tx = (tx_NX == ht_NX) && (tx_NY == ht_NY);
        w->cell_tx = (tx_NX == ht_NX-1) && (tx_NY == ht_NY-1);
        demand(w->node_tx || w->cell_tx, "texture file has wrond size");
        w->texturize = TRUE;
      }
    else
      { w->tx = NULL; 
        w->texturize = w->node_tx = w->cell_tx = FALSE;
      }
    
    /* Nominal radius of data: */
    double SZX = (double)ht_NX;
    double SZY = (double)ht_NY;
    w->RAD = sqrt((SZX*SZX + SZY*SZY)/2);
    
    /* Do not show reference planes initially: */
    w->refPlane = FALSE;
    w->orgPlane = FALSE;
    
    /* Initial observer's position: */
    w->azimuth = 0;
    w->elevation = 30;
    w->distance = (GLfloat)(2*w->RAD);
    
    return w;
  }

fvw_state_range_t fvw_state_range_make(options_range_t *orn, int32_t NC)
  { fvw_state_range_t srn;
    srn.def_vmin = (orn->vauto ? NAN : orn->vmin);
    srn.def_vmax = (orn->vauto ? NAN : orn->vmax);
    srn.vmin = float_vec_new((uint32_t)NC);
    srn.vmax = float_vec_new((uint32_t)NC);
    for (uint32_t c = 0;  c < NC; c++)
      { srn.vmin.e[c] = NAN; 
        srn.vmax.e[c] = NAN;
      }
    return srn;
  }
    
void fvw_get_vrange(uint32_t c, fvw_state_range_t *srn, float_image_t *ht, float *vminP, float *vmaxP)
  { float vmin = srn->vmin.e[c];
    float vmax = srn->vmax.e[c];
    if (isnan(vmin) || isnan(vmax))
      { if (isnan(srn->def_vmin) || isnan(srn->def_vmax))
          { vmin = +INF;
            vmax = -INF;
            float_image_update_sample_range(ht, (int32_t)c, &vmin, &vmax);
            if (vmin > vmax)
              { /* There are no valid samples: */
                vmin = 0; vmax = 0;
              }
          }
        else
          { vmin = srn->def_vmin;
            vmax = srn->def_vmax;
          }
        srn->vmin.e[c] = vmin; 
        srn->vmax.e[c] = vmax; 
      }
    assert(! isnan(vmin));
    assert(! isnan(vmax));
    (*vminP) = vmin;
    (*vmaxP) = vmax;
  }

float_image_t *fvw_read_image(char *name, bool_t isMask, int32_t *NCP, int32_t *NXP, int32_t *NYP)
  { float_image_t *fim = NULL;
    /* Decide type by file name: */
    uint32_t len = (uint32_t)strlen(name);
    if ((strcmp(name,"-") == 0) || ((len > 4) && (strcmp(name+len-4,".fni") == 0)))
      { /* Float image file: */
        fim = fvw_read_float_image(name);
      }
    else if ((len > 4) && (name[len-4] == '.') && (name[len-3] == 'p') && (name[len-1] == 'm'))
      { /* PBM/PGM/PPM file: */
        fim = float_image_read_pnm_named(name, isMask, 1.0, 0.0, TRUE, TRUE, TRUE);
      }
    else
      { fprintf(stderr, "file name = \"%s\"\n", name);
        demand(FALSE, "invalid file name");
      }
    (*NCP) = (int32_t)fim->sz[0];
    (*NXP) = (int32_t)fim->sz[1];
    (*NYP) = (int32_t)fim->sz[2];
    demand (((*NCP) > 0) && ((*NXP) > 0) && ((*NYP) > 0), "image is empty\n");
    return fim;
  }

float_image_t *fvw_read_float_image(char *name)
  { FILE *rd = open_read(name, TRUE);
    float_image_t *fim = float_image_read(rd);
    if (rd != stdin) { fclose(rd); }
    return fim;
  }

void fvw_paint_everything(fvw_state_t *w)
  {
    /* Set the GL graphics modes for depth painting: */
    glShadeModel(GL_SMOOTH);
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    
    /* Get  the height image dimensons: */
    int32_t ht_NC, ht_NX, ht_NY;
    float_image_get_size(w->ht, &ht_NC, &ht_NX, &ht_NY);
    
    /* Get the channel {c} to display: */
    uint32_t c = w->channel;
    
    /* Recompute the height scale range if needed: */
    float ht_vmin, ht_vmax;
    fvw_get_vrange(c, &(w->ht_range), w->ht, &(ht_vmin), &(ht_vmax));
    double ht_vdif = ((double)ht_vmax) - ((double)ht_vmin);
    double ht_vmid = (((double)ht_vmin) + ((double)ht_vmax))/2;
    
    /* Recover the extra height scale factor for this channel: */
    double ht_mag = w->ht_mag.e[c];
    if (isnan(ht_mag))
      { ht_mag = w->def_ht_mag;
        w->ht_mag.e[c] = ht_mag;
      }
    assert (! isnan(ht_mag));
    
    /* Make {ht_mag} the default height scale factor for any other "new" channels: */
    w->def_ht_mag = ht_mag;
    
    /* Compute the total height scale factor: */
    double ht_scale = (ht_vdif == 0.0 ? 1.0 : ht_mag/ht_vdif);
    
    /* Recompute the color scale range if needed: */
    float cm_vmin, cm_vmax;
    fvw_get_vrange(c, &(w->cm_range), w->ht, &(cm_vmin), &(cm_vmax));
    
    /* Compute the center of attention: */
    GLfloat ctrx = (GLfloat)((double)ht_NX/2);
    GLfloat ctry = (GLfloat)((double)ht_NY/2);
    GLfloat ctrz = (GLfloat)(ht_scale * ht_vmid);

    /* Compute observer's position: */
    double az =  M_PI*w->azimuth/180;   /* Observer's azimuth in radians. */
    double ev =  M_PI*w->elevation/180; /* Observer's elevation in radians. */
    double ca = cos(az), sa = sin(az);  
    double ce = cos(ev), se = sin(ev);
    double R = w->distance;

    GLfloat obsx = (GLfloat)(ctrx + R * ca * ce);
    GLfloat obsy = (GLfloat)(ctry + R * sa * ce);
    GLfloat obsz = (GLfloat)(ctrz + R * se);

    glPushMatrix();
    
    /* Set the view transformation matrix: */
    gluLookAt( obsx, obsy, obsz, ctrx, ctry, ctrz, 0.0, 0.0, 1.0);

    /* Place the lights: */
    fvw_set_lights(w);
    
    /* Paint the reference planes, if any: */
    if (w->refPlane)
      { double eps = fmax(0.00001*w->RAD, 0.001*ht_vdif);
        float zmin = (float)(ht_scale*ht_vmin - eps);
        fvw_paint_reference_plane(w, zmin);
        float zmax = (float)(ht_scale*ht_vmax + eps);
        fvw_paint_reference_plane(w, zmax);
      }
    if (w->orgPlane)
      { float zorg = 0;
        fvw_paint_reference_plane(w, zorg);
      }

    /* Paint the terrain: */
    fvw_paint_height_map
      ( w->ht, c, ht_scale, cm_vmin, cm_vmax, w->hist,
        (w->texturize ? w->tx : NULL), w->node_tx, w->cell_tx
      );
    
    glPopMatrix();
    
    glDisable(GL_DEPTH_TEST);
  }    

void fvw_set_perspective(fvw_state_t *w)
  {
    /* Set viewport: */
    glViewport(0, 0, (GLint)w->window_HSize, (GLint)w->window_VSize);

    /* Set perspective matrix: */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float angheight = 40.0; /* Angular height of window (degrees) */
    float aspect = (float)w->window_HSize/(float)w->window_VSize;
    gluPerspective(angheight, aspect, 0.1, w->distance + 2.0*w->RAD);
    glMatrixMode(GL_MODELVIEW);
  }

void fvw_set_lights(fvw_state_t *w)
  {
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);

    double ambi = 0.1;                /* Amount of ambient light. */
    double spec = 0.0;                /* Amount of shine-scatterable light. */
    double diff = 1.0 - ambi - spec;  /* Amount of diffusible light. */

    if (ambi > 0)
      { GLfloat light_ambient[] = { (GLfloat)ambi, (GLfloat)ambi, (GLfloat)ambi, (GLfloat)1.0 };
        glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
      }

    if (diff > 0)
      { GLfloat light_diffuse[] = { (GLfloat)diff, (GLfloat)diff, (GLfloat)diff, (GLfloat)1.0 };
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
      }

    if (spec > 0)
      { GLfloat light_specular[] = { (GLfloat)spec, (GLfloat)spec, (GLfloat)spec, (GLfloat)1.0 };
        glLightfv(GL_LIGHT0, GL_SPECULAR,  light_specular);
      }

    /* Choose light's position relative to observer: */
    double rev = w->elevation/90.0;  /* Observer's rel elevation, in {[-1_+1]}. */
    double daz = +30;     /* Azimuth rel to observer (degrees). */
    double dev = rev*30;  /* Elevation rel to observer (degrees). */
    /* Compute light  direction vector: */
    double az =  M_PI*(w->azimuth + daz)/180;   /* Light's azimuth (radians). */
    double ev =  M_PI*(w->elevation + dev)/180; /* Observer's elevation (radians). */
    double ca = cos(az), sa = sin(az);  
    double ce = cos(ev), se = sin(ev);
    GLfloat dirx = (GLfloat)(ca * ce);
    GLfloat diry = (GLfloat)(sa * ce);
    GLfloat dirz = (GLfloat)se;
    
    /* Place light at infinity: */
    GLfloat pos[] = { dirx, diry, dirz, 0.0 };
    glLightfv(GL_LIGHT0, GL_POSITION, pos);
  }

void fvw_paint_height_map
  ( float_image_t *ht, 
    uint32_t c, 
    double ht_scale, 
    float cm_vmin,
    float cm_vmax,
    bool_t hist,
    float_image_t *tx,
    bool_t node_tx,
    bool_t cell_tx
  )
  {
    if (ht == NULL) { return; }
    if (fvw_debug_paint) { fprintf(stderr, "+ %s\n", __FUNCTION__); }

    /* Get  the height image dimensons: */
    int32_t ht_NC, ht_NX, ht_NY;
    float_image_get_size(ht, &ht_NC, &ht_NX, &ht_NY);

    /* Set surface finish: */
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    
    if (hist)
      { /* Histogram style, paint with height-derived colors, no skirt: */
        fvw_paint_self_colored_hist_height_map(ht, c, ht_scale, cm_vmin, cm_vmax, FALSE);
      }
    else if (tx == NULL)
      { /* No texture map, paint with height-derived colors: */
        fvw_paint_self_colored_height_map(ht, c, ht_scale, cm_vmin, cm_vmax);
      }
    else if (cell_tx)
      { /* Texture colors are associated with grid cells: */
        fvw_paint_cell_colored_height_map(ht, c, ht_scale, tx);
      }
    else if (node_tx)
      { /* Texture colors are associated with grid corners: */ 
        fvw_paint_node_colored_height_map(ht, c, ht_scale, tx);
      }
    else
      { assert(FALSE); }
    glDisable(GL_COLOR_MATERIAL);
    if (fvw_debug_paint) { fprintf(stderr, "- %s\n", __FUNCTION__); }
  }

void fvw_paint_reference_plane(fvw_state_t *w, float z)
  {
    /* Get the height image dimensons: */
    int32_t ht_NC, ht_NX, ht_NY;
    float_image_get_size(w->ht, &ht_NC, &ht_NX, &ht_NY);
    
    /* X and Y extents are those of the bitmap: */
    float xmin = 0.0;
    float xmax = (float)ht_NX;
    float ymin = 0.0;
    float ymax = (float)ht_NY;
    
    glColorMaterial(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);

    glColor3f(0.650f, 0.975f, 0.925f);
    
    glBegin(GL_QUADS);
    glVertex3f(xmin, ymin, z);
    glVertex3f(xmax, ymin, z);
    glVertex3f(xmax, ymax, z);
    glVertex3f(xmin, ymax, z);
    glEnd();
    
    glDisable(GL_COLOR_MATERIAL);
  }

/* GL EVENT-HANDLING METHODS */

void fvw_display_method(void)
  {
    if (fvw_debug_GL) { fprintf(stderr, "+ %s\n", __FUNCTION__); }
    fvw_state_t *w = fvw_state;

    /* Clear the window: */
    glClearColor(0.750, 0.750, 0.750, 1.000); /* Light gray, opaque. */
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

    /* Paint the terrain and stuff: */
    fvw_paint_everything(w);
    
    /* Display the result: */
    glutSwapBuffers();
    if (fvw_debug_GL) { fprintf(stderr, "- %s\n", __FUNCTION__); }
  }

void fvw_reshape_method(int32_t width, int32_t height)
  {
    if (fvw_debug_GL) { fprintf(stderr, "+ %s\n", __FUNCTION__); }

    fvw_state_t *w = fvw_state;
    
    /* Save window size: */
    w->window_HSize = width;
    w->window_VSize = height;
    
    /* Adjust the perspective parameters to the window's current size: */
    fvw_set_perspective(w);

    if (fvw_debug_GL) { fprintf(stderr, "- %s\n", __FUNCTION__); }
  }

void fvw_keyboard_method(unsigned char key, int32_t x, int32_t y)
  {
    if (fvw_debug_GL) { fprintf(stderr, "+ %s\n", __FUNCTION__); }

    fvw_state_t *w = fvw_state;
    
    switch (key) 
      { 
        case 27: 
        case 'q':
        case 'Q':
          /* Quit: */
          exit(0);
          break;

        case 't':
          /* Toggle texturization: */
          w->texturize = !w->texturize;
          glutPostRedisplay();
          break;

        case 'p':
          /* Toggle reference plane: */
          w->refPlane = !w->refPlane;
          glutPostRedisplay();
          break;

        case 'o':
          /* Toggle {Z=0} plane: */
          w->orgPlane = !w->orgPlane;
          glutPostRedisplay();
          break;

        case 'h':
          /* Toggle histogram mode: */
          w->hist = !w->hist;
          glutPostRedisplay();
          break;

        case 'Z':
          /* Zoom in, down to a limit: */
          w->distance = (GLfloat)fmax(w->distance*fvw_zoom_step_ratio, w->RAD/16);
          glutPostRedisplay();
          break;

        case 'z':
          /* Zoom out, up to a limit: */
          w->distance = (GLfloat)fmin(w->distance/fvw_zoom_step_ratio, w->RAD*16);
          glutPostRedisplay();
          break;

        case 'S':
          /* Increase extra height scale factor: */
          w->ht_mag.e[w->channel] *= fvw_ht_scale_step_ratio;
          glutPostRedisplay();
          break;

        case 's':
          /* Reduce extra height scale factor: */
          w->ht_mag.e[w->channel] /= fvw_ht_scale_step_ratio;
          glutPostRedisplay();
          break;

        case 'c':
          /* Cycle through channels of height map: */
          { uint32_t NC = (uint32_t)w->ht->sz[0]; 
            w->channel = (w->channel + 1) % NC;
          }
          glutPostRedisplay();
          break;

        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
          /* Select named channel of height map: */
          { int32_t c = key - '0';
            assert((c >= 0) && (c <= 9));
            int32_t NC = (int32_t)w->ht->sz[0];
            w->channel = (uint32_t)(c < NC ? c : NC-1);
          }
          glutPostRedisplay();
          break;
      }

    if (fvw_debug_GL) { fprintf(stderr, "- %s\n", __FUNCTION__); }
  }

void fvw_passivemouse_method( int32_t x, int32_t y)
  {
    fvw_state_t *w = fvw_state;
    
    w->mouse_x = x;
    w->mouse_y = y;
  }

void fvw_activemouse_method(int32_t x, int32_t y)
  {
    fvw_state_t *w = fvw_state;
    
    int32_t xvec = w->mouse_x - x;
    int32_t yvec = w->mouse_y - y;

    double az = (double)w->azimuth + xvec;
    double el = (double)w->elevation - yvec;
    el = fmax(-88.0, fmin(+88.0, el));
    
    w->elevation = (GLfloat)el;
    w->azimuth = (GLfloat)az;

    w->mouse_x = x;
    w->mouse_y = y;

    glutPostRedisplay();
  }

void fvw_special_method(int32_t key, int32_t x, int32_t y)
  {
    if (fvw_debug_GL) { fprintf(stderr, "+ %s\n", __FUNCTION__); }

    fvw_state_t *w = fvw_state;
    
    switch (key)
      { 
        case GLUT_KEY_UP:
          w->elevation = (GLfloat)fmin(w->elevation + 2, +88.0);
          glutPostRedisplay();
          break;

        case GLUT_KEY_DOWN:
          w->elevation = (GLfloat)fmax(w->elevation - 2, -88.0);
          glutPostRedisplay();
          break;

        case GLUT_KEY_LEFT:
          w->azimuth = w->azimuth + 2;
          glutPostRedisplay();
          break;

        case GLUT_KEY_RIGHT:
          w->azimuth = w->azimuth - 2;
          glutPostRedisplay();
          break;
      }

    if (fvw_debug_GL) { fprintf(stderr, "- %s\n", __FUNCTION__); }
  }

/* COMMAND LINE PARSING */

options_t *fvw_parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
    
    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    /* Parse the texture map name: */
    if (argparser_keyword_present(pp, "-texture"))
      { o->texture_file = argparser_get_next(pp); }
    else
      { o->texture_file = NULL; }

    /* Parse the initial channel index to display: */
    if (argparser_keyword_present(pp, "-channel"))
      { o->channel = (uint32_t)argparser_get_next_int(pp, 0, float_image_max_size-1); }
    else
      { o->channel = 0; }

    /* Parse the reference scaling range: */
    o->range = fvw_parse_options_range(pp, "-range");
    
    /* Parse the initial scale to display: */
    if (argparser_keyword_present(pp, "-scale"))
      { o->scale = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); }
    else
      { o->scale = 1.00; }

    /* Parse the integer decoding options: */
    if (argparser_keyword_present(pp, "-isMask"))
      { o->isMask = argparser_get_next_bool(pp); }
    else
      { o->isMask = FALSE; }
    
    /* Parse the histogram style option: */
    if (argparser_keyword_present(pp, "-hist"))
      { o->hist = argparser_get_next_bool(pp); }
    else
      { o->hist = FALSE; }
    
    /* Parse the colormapping range: */
    o->colorize = fvw_parse_options_range(pp, "-colorize");

    /* Go to the positional arguments: */
    argparser_skip_parsed(pp);
    
    /* Parse the height map name and its selected channel: */
    if ((argparser_next(pp) != NULL) && (! argparser_next_is_keyword(pp)))
      { o->height_file = argparser_get_next(pp); }
    else
      { o->height_file = "-"; }
    
    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

options_range_t fvw_parse_options_range(argparser_t *pp, char *keyword)
  { options_range_t orn;
    orn.vmin = NAN; /* Default. */
    orn.vmax = NAN; /* Default. */
    if (argparser_keyword_present(pp, keyword))
      { if (argparser_keyword_present_next(pp, "auto"))
          { orn.vauto = TRUE; }
        else
          { orn.vauto = FALSE;
            orn.vmin = (float)argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
            orn.vmax = (float)argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
          }
      }
    else
      { orn.vauto = TRUE; }
    return orn;
  }
