#define PROG_NAME "fni_view"
#define PROG_DESC "3D interactive visualization of a height map"
#define PROG_VERS "2.0"

#define PROG_C_COPYRIGHT "Copyright © 2005 Universidade Estadual Fluminense (UFF)."

/* Last edited on 2025-01-23 14:21:29 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -range {VMIN} {VMAX} | -range auto ] \\\n" \
  "    [ -scale {VSCALE} ] \\\n" \
  "    [ -isMask {ISMASK} ] \\\n" \
  "    [ -colorize {VMIN} {VMAX} | -colorize auto ] \\\n" \
  "    [ -channel {CHAN} ] \\\n" \
  "    [ -txFile {TXFILE} | -hist {HISTFLAG}  ] \\\n" \
  "    [ -txChannels {TXCHAN}... ] \\\n" \
  "    [ -verbose ] \\\n" \
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
  "COLORIZATION\n" \
  "  The terrain or histogram is colorized either by otaining, for each height" \
  " value, either an RGB texture" \
  " value, or a gray texture value and converting it" \
  " to an RGB color by an internal colorizing" \
  " palette.   The texture values (RGB or gray) can be" \
  " obtained from specific channels of" \
  " a separate texture image or from the" \
  " currently displayed channel of the input image.  See" \
  " the \"-txFile\" and \"-txChannels\" options below.\n" \
  "\n" \
  "  If a separate texture image is specified, it must have the same" \
  " size as the input image, image, or (in terrain mode only) exactly one column" \
  " and one row less than it.  In the" \
  " first case, the texture pixel in column {x}, row {y} is taken" \
  " as the color of the terrain above the unit square with" \
  " diagonal {(x-1/2,y-1/2)} and {(x+1/2,y+1/2)}.  In the" \
  " second case, the pixel in column {x}, row {y} of the texture map is" \
  " taken as the color of the unit-side square whose diagonal" \
  " corners are {(x,y)} and {(x+1,y+1)}.\n" \
  "\n" \
  "  The input file (and the texture file, if any) may" \
  " be in the PNM (PBM/PGM/PPM) format" \
  " (see {uint16_image.h}) or in the float-valued multi-channel" \
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
  "  -txFile {TXFILE}\n" \
  "    Specifies the name of the FNI file containing the texture" \
  " map.  See also the interactive keyboard command 't'.  If {TXFILE}" \
  " is \"-\", reads the texture map from standard input too (after the height" \
  " map, if applicable).  If omitted, the input image itself is used as the texture image.\n" \
  "\n" \
  "  -txChannels {TX_CHAN}\n" \
  "  -txChannels {TX_CHAN0} {TX_CHAN1} {TX_CHAN2}\n" \
  "    This optional parameter specifies either one or three channels of the" \
  " texture file (or of the input file, if no \"-txFile\" was specified).  If only" \
  " one channel is specified, its sample values are converted to RGB colors using" \
  " an internal colorizing palette.  If \"-txFile\" is given but \"-txChannels\" is" \
  " omitted, assumes \"-txChannels 0 1 2\" if the texture file has three or more" \
  " channels, and \"-txChannels 0\" if it has fewer than three.  If" \
  " both \"-txFile\" and \"-txChannels\" are omitted, uses the channel" \
  " of the input image that is currently being displayed.\n" \
  "\n" \
  "  -hist {HISTFLAG}\n" \
  "    If {HISTFLAG} is true (\"T\" or 1), the map is to be drawn" \
  " as an histogram, with each pixel as a horizontal square and vertical" \
  " walls between adjacent pixels. In this case the \"-txFile\" option" \
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
  " omitted or turned off, or has a single channel.  In the first option, samples" \
  " in the range {[VMIN _ VMAX]} are mapped to colors of an appropriate" \
  " pseudocolor scale; samples outside that range are mapped to the end" \
  " points of the scale.  The \"auto\" option will set {VMIN} and {VMAX} to" \
  " the maximum and minimum" \
  " samples in the current channel.  The default is \"-colorize auto\".  The" \
  " color mapping is not affected by the \"-range\" and \"-scale\" arguments" \
  " or the 's' and 'S' command keys.\n" \
  "\n" \
  "  -verbose\n" \
  "    If present, the program writes details of the file input and processing.\n" \
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
  "  jan/2025 by Jorge Stolfi: Added \"-txChannels\" option.\n" \
  "  jan/2025 by Jorge Stolfi: Added \"-verbose\" option.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <values.h>
#include <assert.h>
#include <string.h>

#define __APPLE__ 0
#include <GL/glu.h>
#include <GL/glut.h>

#include <argparser.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_read_pnm.h>
#include <sample_conv.h>
#include <frgb.h>
#include <bool.h>
#include <jsfile.h>
#include <r3.h>

#include <fvw_GL.h>
#include <fvw_paint.h>
#include <fvw_paint_height_map.h>

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

typedef struct options_t
  { char* height_file;        /* Name of height map. */
    /* Height scaling options: */
    options_range_t range;    /* Height range options for scaling. */
    double scale;             /* Additional height scale factor. */
    /* Colorizing options: */
    char* txFile;             /* Name of texture map, or NULL if not specified. */
    int32_t txChannels[3];    /* Texture map channels to use, or {-1} if not specified. */
    options_range_t colorize; /* Height range options for color mapping. */
    /* Other options: */
    uint32_t channel;         /* Initial channel to display. */
    bool_t isMask;            /* True if input samples are to be converted as in a mask. */
    bool_t hist;              /* True specifies histogram-like plot. */
    bool_t verbose;           /* True asks for diagnostics of files etc. */
  } options_t;

/* GLOBAL VARIABLES USED BY GL METHODS */

typedef struct fvw_state_range_t
  { float def_vmin;        /* Default {vmin} for new channels or {NAN} if auto. */
    float def_vmax;        /* Default {vmax} for new channels or {NAN} if auto. */
    float_vec_t vmin;      /* High pix values per channel or {NAN} if undef (mutable). */
    float_vec_t vmax;      /* Low pix values per channel or {NAN} if undef (mutable). */
  } fvw_state_range_t;
  /* Specifies range parameters for height scaling or colormapping. */

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
    fvw_texture_t *tx;   /* Texture information. */
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
    
    /* Verbosity and diagnostics: */
    bool_t verbose; /* Print informative messages. */
  } fvw_state_t;

static fvw_state_t *fvw_state = NULL; /* The state displayed in the GL window. */

/* INTERNAL PROTOTYPES */
 
options_range_t fvw_parse_options_range(argparser_t *pp, char *keyword);

fvw_state_range_t fvw_state_range_make(options_range_t *orn, int32_t NC);
  /* Creates a {fvw_state_range_t} from the user options {orn} for
    a height map image with {NC} channels. */

fvw_state_t *fvw_create_state(options_t *o);
  /* Creates and initializes a window data record, suitable for use
    used by the methods {fvw_display_method} etc. */

fvw_texture_t* fvw_make_texture_record(options_t *o, float_image_t *ht);
  /* Creates a {fvw_texture_t} from options including {o.txFile}, {o.txChannels}, {o.colorize}.
    If a texture file was specified, checks if its size is compatible with the height map {ht}.
    If no texture file was specified, but texture channels were, uses the height map {ht}.
    If neither texture file nor texture channels were specified, returns {NULL}. */
  
void fvw_get_vrange(uint32_t c, fvw_state_range_t *srn, float_image_t *ht, float *vminP, float *vmaxP);  
  /* Obtains the height scaling or colormapping pixel range {*vminP} and {*vmaxP} for channel {c}
    of the image {ht}.  If the range has been previously computed and saved in the record {srn},
    takes it from there.  Otherwise, if {srn} has a user-specified range, uses that range.
    Otherwise recomputes the range from the pixel values of {ht}. Either way, saves these
    values in the {srn} record for future use. */

void  fvw_get_sample_range_from_image(float_image_t *timg, int32_t c, float *vmin_P, float *vmax_P);
  /* Gets the min and max sample values in channel {c} of {timg}, fudging
    it if is trivial or empty. */

float_image_t *fvw_read_image(char *name, bool_t isMask, int32_t *NCP, int32_t *NXP, int32_t *NYP, bool_t verbose);
  /* Reads a float image, in the format of {float_image_write} (if the
    {name} ends in ".fni") or in PBM/PGM/PPM format (if {name} ends in
    ".pbm", ".pgm", ".ppm", ".pnm"). If {name} is "-", reads from
    {stdin} assuming {float_image_write} format. The procedure also
    sets {*NCP,*NXP,*NYP} to the number of channels, columns, and rows.*/

float_image_t *fvw_read_float_image(char *name, bool_t verbose);
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
    int32_t NC_ht, NX_ht, NY_ht;
    w->ht = fvw_read_image(o->height_file, o->isMask, &NC_ht, &NX_ht, &NY_ht, o->verbose);
    if (o->verbose) { fprintf(stderr, "height map has %d channels, %d columns, %d rows\n", NC_ht, NX_ht, NY_ht); }
    demand(NC_ht > 0, "invalid num of channels in height map");
    
    /* The initial channel is user-specified: */
    w->channel = o->channel;
    if ((w->channel < 0) || (w->channel >= NC_ht))
      { fprintf(stderr, "invalid height map channel %d, using 0\n", w->channel);
        w->channel = 0;
      }

    /* Pixel ranges for height scaling: */
    w->ht_range = fvw_state_range_make(&(o->range), NC_ht);

    /* Extra height scaling factors: */
    w->def_ht_mag = o->scale;
    w->ht_mag = double_vec_new((uint32_t)NC_ht);
    for (uint32_t c = 0;  c < NC_ht; c++) { w->ht_mag.e[c] = NAN; }

    /* Pixel ranges for color mapping, if not texturized: */
    w->cm_range = fvw_state_range_make(&(o->colorize), NC_ht);
    
    /* The histogram-style flag is initially obtained from the command line: */
    w->hist = o->hist;
    w->tx = fvw_make_texture_record(o, w->ht);
    w->texturize = (w->tx != NULL); /* If there is texture, start texturized. */
    
    /* Nominal radius of data: */
    double SZX = (double)NX_ht;
    double SZY = (double)NY_ht;
    w->RAD = sqrt((SZX*SZX + SZY*SZY)/2);
    
    /* Do not show reference planes initially: */
    w->refPlane = FALSE;
    w->orgPlane = FALSE;
    
    /* Initial observer's position: */
    w->azimuth = 0;
    w->elevation = 30;
    w->distance = (GLfloat)(2*w->RAD);
    
    w->verbose = o->verbose;
    
    return w;
  }
    
fvw_texture_t* fvw_make_texture_record(options_t *o, float_image_t *ht)
  { 
    /* Get size of height map: */
    int32_t NC_ht, NX_ht, NY_ht;
    float_image_get_size(ht, &NC_ht, &NX_ht, &NY_ht);

    fvw_texture_t *tx = NULL;
    int32_t NC_tx, NX_tx, NY_tx; /* Channels, cols and rows of texture map. */
    
    if (o->txFile != NULL)
      { /* A texture file was specified: */
        tx = talloc(1, fvw_texture_t);
        /* Read texture map, if any,and check its size: */
        bool_t tx_isMask = FALSE; /* Assume that texture has smooth pixel distr. */
        tx->timg = fvw_read_image(o->txFile, tx_isMask, &NC_tx, &NX_tx, &NY_tx, o->verbose);
        if (o->verbose) { fprintf(stderr, "texture map has %d channels, %d columns, %d rows\n", NC_tx, NX_tx, NY_tx); }
        if ((NX_tx == NX_ht) && (NY_tx == NY_ht))
          { tx->cell = FALSE; }
        else if ((NX_tx == NX_ht-1) && (NY_tx == NY_ht-1))
          { tx->cell = TRUE; }
        else
          { demand(FALSE, "texture file has wrong size"); }
      }
    else if (o->txChannels[0] != -1)
      { /* Texture channels were specified but no texture file was, use height map: */
        if (o->verbose) { fprintf(stderr, "using the height map as texture image\n"); }
        tx = talloc(1, fvw_texture_t);
        tx->timg = ht; 
        NC_tx = NC_ht; NX_tx = NX_ht; NY_tx = NY_ht;
        tx->cell = FALSE;
      }
    else
      { /* Neither a texture file not texture channels were specified, no texture: */
        tx = NULL; 
        NC_tx = -1; NX_tx = -1; NY_tx = -1;
      }
      
    if (tx != NULL)
      { /* Set the texture channels to use: */
        for (int32_t k = 0; k < 3; k++)
          { int32_t ch = o->txChannels[k];
            demand((ch == -1) || ((ch >= 0) && (ch < NC_tx)), "invalid texture channel");
            tx->chans[k] = ch;
          }

        if (tx->chans[0] == -1)
          { /* There is a tx image but channels were not specified, provide defaults: */
            tx->chans[0] = 0;
            if (NC_tx >= 3) { tx->chans[1] = 1; tx->chans[2] = 2; }
          }
          
        demand(tx->chans[0] != -1, "must specify at least 1 texture channel");
        demand((tx->chans[1] == -1) == (tx->chans[2] == -1), "must specify 1 or 3 texture channels");
        
        if (tx->chans[1] == -1)
          { /* Define {tx.vmin,tx.vmax} for pseudocolor generation: */
            if (o->colorize.vauto)
              { fvw_get_sample_range_from_image(tx->timg, tx->chans[0], &(tx->vmin), &(tx->vmax)); }
            else
              { tx->vmin = o->colorize.vmin;
                tx->vmax = o->colorize.vmax;
              }
            if (o->verbose) { fprintf(stderr, "texture map colorization range = [ %15.7e _ %15.7e ]\n", tx->vmin, tx->vmax); }
          }
        else
          { /* NO need for {tx.vmin,tx.vmax} */ }
      }
    return tx;
  }
    
void fvw_get_vrange(uint32_t c, fvw_state_range_t *srn, float_image_t *ht, float *vminP, float *vmaxP)
  { float vmin = srn->vmin.e[c];
    float vmax = srn->vmax.e[c];
    if (isnan(vmin) || isnan(vmax))
      { if (isnan(srn->def_vmin) || isnan(srn->def_vmax))
          { fvw_get_sample_range_from_image(ht, (int32_t)c, &vmin, &vmax); }
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
              
void  fvw_get_sample_range_from_image(float_image_t *timg, int32_t c, float *vmin_P, float *vmax_P)
  { 
    float vmin = +INF;
    float vmax = -INF;
    float_image_update_sample_range(timg, c, &(vmin), &(vmax));
    if (vmin > vmax)
      { vmin = -1.0; vmax = +1.0; }
    else 
      { float vdel = vmax - vmin;
        float vmag = fmaxf(fabsf(vmax), fabsf(vmin));
        float eps = (float)(fmax(1.0e-5*vmag, 1.0e-5));
        if (vdel < 2*eps) { vmin -= eps; vmax += eps; }
      }
    (*vmin_P) = vmin;
    (*vmax_P) = vmax;
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

float_image_t *fvw_read_image(char *name, bool_t isMask, int32_t *NCP, int32_t *NXP, int32_t *NYP, bool_t verbose)
  { float_image_t *fim = NULL;
    /* Decide type by file name: */
    uint32_t len = (uint32_t)strlen(name);
    if ((strcmp(name,"-") == 0) || ((len > 4) && (strcmp(name+len-4,".fni") == 0)))
      { /* Float image file: */
        fim = fvw_read_float_image(name, verbose);
      }
    else if ((len > 4) && (name[len-4] == '.') && (name[len-3] == 'p') && (name[len-1] == 'm'))
      { /* PBM/PGM/PPM file: */
        fim = float_image_read_pnm_named(name, isMask, 1.0, 0.0, TRUE, verbose, verbose);
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

float_image_t *fvw_read_float_image(char *name, bool_t verbose)
  { FILE *rd = open_read(name, verbose);
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
    int32_t NC_ht, NX_ht, NY_ht;
    float_image_get_size(w->ht, &NC_ht, &NX_ht, &NY_ht);
    
    /* Get the channel {c} to display: */
    uint32_t c = w->channel;
    
    /* Recompute the height scale range if needed: */
    float ht_vmin, ht_vmax;
    fvw_get_vrange(c, &(w->ht_range), w->ht, &(ht_vmin), &(ht_vmax));
    if (w->verbose) { fprintf(stderr, "channel range = [ %15.7e _ %15.7e ]\n", ht_vmin, ht_vmax); }
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
    if (w->verbose) { fprintf(stderr, "height scale factor = %15.7e\n", ht_scale); }
    
    /* Recompute the color scale range if needed: */
    float cm_vmin, cm_vmax;
    fvw_get_sample_range_from_image(w->ht, (int32_t)c, &(cm_vmin), &(cm_vmax));
    
    /* Compute the center of attention: */
    GLfloat ctrx = (GLfloat)((double)NX_ht/2);
    GLfloat ctry = (GLfloat)((double)NY_ht/2);
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
      ( w->ht, c, ht_scale,
        w->hist, 
        (w->texturize ? w->tx : NULL),
        cm_vmin, cm_vmax
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


void fvw_paint_reference_plane(fvw_state_t *w, float z)
  {
    /* Get the height image dimensons: */
    int32_t NC_ht, NX_ht, NY_ht;
    float_image_get_size(w->ht, &NC_ht, &NX_ht, &NY_ht);
    
    /* X and Y extents are those of the bitmap: */
    float xmin = 0.0;
    float xmax = (float)NX_ht;
    float ymin = 0.0;
    float ymax = (float)NY_ht;
    
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
          if (w->verbose) { fprintf(stderr, "viewing distance = %15.7e\n", w->distance); }
          glutPostRedisplay();
          break;

        case 'z':
          /* Zoom out, up to a limit: */
          w->distance = (GLfloat)fmin(w->distance/fvw_zoom_step_ratio, w->RAD*16);
          if (w->verbose) { fprintf(stderr, "viewing distance = %15.7e\n", w->distance); }
          glutPostRedisplay();
          break;

        case 'S':
          /* Increase extra height scale factor: */
          w->ht_mag.e[w->channel] *= fvw_ht_scale_step_ratio;
          if (w->verbose) { fprintf(stderr, "extra scale factor = %15.7e\n", w->ht_mag.e[w->channel]); }
          glutPostRedisplay();
          break;

        case 's':
          /* Reduce extra height scale factor: */
          w->ht_mag.e[w->channel] /= fvw_ht_scale_step_ratio;
          if (w->verbose) { fprintf(stderr, "extra scale factor = %15.7e\n", w->ht_mag.e[w->channel]); }
          glutPostRedisplay();
          break;

        case 'c':
          /* Cycle through channels of height map: */
          { uint32_t NC = (uint32_t)w->ht->sz[0]; 
            w->channel = (w->channel + 1) % NC;
            if (w->verbose) { fprintf(stderr, "showing channel %d\n", w->channel); }
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
            if (w->verbose) { fprintf(stderr, "showing channel %d\n", w->channel); }
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
    
    options_t *o = (options_t *)malloc(sizeof(options_t)); 
    
    if (argparser_keyword_present(pp, "-txFile"))
      { o->txFile = argparser_get_next(pp); }
    else
      { o->txFile = NULL; }
    
    for (int32_t k = 0; k < 3; k++) { o->txChannels[k] = -1; }
    if (argparser_keyword_present(pp, "-txChannels"))
      { o->txChannels[0] = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX);
        if (argparser_next_is_number(pp))
          { o->txChannels[1] = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX);
            o->txChannels[2] = (int32_t)argparser_get_next_int(pp, 0, INT32_MAX);
          }
      }

    if (argparser_keyword_present(pp, "-channel"))
      { o->channel = (uint32_t)argparser_get_next_int(pp, 0, float_image_max_size-1); }
    else
      { o->channel = 0; }

    o->range = fvw_parse_options_range(pp, "-range");
    
    if (argparser_keyword_present(pp, "-scale"))
      { o->scale = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX); }
    else
      { o->scale = 1.00; }

    if (argparser_keyword_present(pp, "-isMask"))
      { o->isMask = argparser_get_next_bool(pp); }
    else
      { o->isMask = FALSE; }
    
    if (argparser_keyword_present(pp, "-hist"))
      { o->hist = argparser_get_next_bool(pp); }
    else
      { o->hist = FALSE; }
    
    o->colorize = fvw_parse_options_range(pp, "-colorize");
    
    o->verbose = argparser_keyword_present(pp, "-verbose");

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
