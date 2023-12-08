#define PROG_NAME "ppmoquant"
#define PROG_DESC "quantize two PPM files into a partially transparent PPM"
#define PROG_VERS "1.0"

/* Last edited on 2021-07-17 23:50:21 by jstolfi */

#define ppmoquant_C_COPYRIGHT \
  "Copyright © 1996 by the State University of Campinas (UNICAMP) - Rewritten version"  

/* TO DO: !!! Currently has syntax errors -- Update calls and fix inconsistencies. */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "    [ -fs | -floyd | -nofs | -nofloyd ] \\\n" \
  "    -transparent {RVAL} {GVAL} {BVAL} \\\n" \
  "    -map {MAP_PPMFILE} \\\n" \
  "    [ -verbose ] \\\n" \
  "    {IN_PPMFILE_0} {IN_PPMFILE_1} \\\n" \
  "    > {OUT_PPMFILE}"

#define PROG_INFO_DESC \
  "  This program is similar to Jef Poskanzer's {ppmquant.c}, but" \
  " takes /two/ images of an object, {IN_PPMFILE_0} and {IN_PPMFILE_1}, which" \
  " are assumed to have been photographed or rendered over black and white" \
  " backgrounds, respectively.  The program then uses a given set of" \
  " colors plus a single `transparent' color to produce a single" \
  " dithered image {OUT_PPMFILE} that best reproduces those two images.\n" \
  "\n"\
  "  The dithered image is always written to {stdout}, and has" \
  " the same {maxval} as the input images or {maxval = 255}, whichever" \
  " is larger.  In regions" \
  " where {IN_PPMFILE_0} is black and {IN_PPMFILE_1} is white, the" \
  " output image will have transparent pixels.  In regions where" \
  " the two input images have the same color, the output image" \
  " will have a mix of opaque colors that approximates that same" \
  " color on the average.  In other regions, opaque and" \
  " transparent pixels are mixed to emulate the object's" \
  " color and transaprency.\n" \
  "\n"\
  "  Since this program allows only fully opaque or fully" \
  " transparent pixels it can only emulate alpha-channel-type" \
  " transparency, like that of a thin layer of opaque dust or" \
  " an opaque screen with tiny holes. It cannot imitate the look" \
  " of objects that have different transmission in different color" \
  " bands, such as colored glass or tea.\n" \
  "\n"\
  "  The set of allowed opaque colors is specified by a" \
  " PPM file, the {MAP_PPMFILE}.  For instance, to quantize" \
  " down to the eight corners of the color cube, you might" \
  " provide the following file as the {MAP_PPMFILE}:\n" \
  "\n" \
  "     P3\n" \
  "     8 1\n" \
  "     255\n" \
  "       0   0   0\n" \
  "     255   0   0\n" \
  "       0 255   0\n" \
  "       0   0 255\n" \
  "     255 255   0\n" \
  "     255   0 255\n" \
  "       0 255 255\n" \
  "     255 255 255\n" \
  "\n" \
  "   If you want to quantize one pixmap to use the colors in another" \
  " one, just use the second one as the mapfile.  You don't have to" \
  " reduce it down to only one pixel of each color, just use it as is.\n" \
  "\n" \
  "   Since the PPM format of the output file does not allow" \
  " transparent pixels or alpha values, the transparent pixels are" \
  " represented by pixels with the specified {COLOR}. This color must" \
  " not occur in the {MAP_PPMFILE}."

#define PROG_INFO_OPTS \
  "  -floyd\n"\
  "  -fs\n" \
  "  -nofloyd\n"\
  "  -nofs\n" \
  "    The optional flag \"-fs\" (or \"-floyd\") enables the Floyd-Steinberg" \
  " error diffusion algorithm, which gives vastly better results" \
  " on images where the unmodified quantization has banding or" \
  " other artifacts, especially when going to a small number of" \
  " colors. The flag \"-nofs\" (or \"-nofloyd\") disables it.  The" \
  " default is \"-nofloyd\".\n" \
  "\n" \
  "  -map {MAP_PPMFILE}\n" \
  "    This mandatory argument specifies the name of the PPM file" \
  " that contains the set of colors to use in the dithering.\n" \
  "\n" \
  "  -transparent {RVAL} {GVAL} {BVAL}\n" \
  "    This mandatory argument specifies that transparent pixels" \
  " shall be represented in the output PPM file by the RGB" \
  " value {(RVAL,GVAL,BVAL)}, with decimal integer components.\n" \
  "\n" \
  "  -verbose\n" \
  "    This option, if present, generates some debugging messages to {stderr}." \
  
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
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  ppmquant(1), ppmvquant(1), ppmquantall(1), pnmdepth(1), ppmdither(1), ppm(5).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created in 1996 by J. Stolfi based on {ppmquant.c} by Jef Poskanzer.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "   1989--1991 Original {ppmquant.c} created by Jef Poskanzer.\n" \
  "   1996-05-18 Opacity logic added by J. Stolfi.\n" \
  "   2009-01-06 Rewritten to use {argparser,jspnm}, etc. by J. Stolfi.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " ppmoquant_C_COPYRIGHT ".\n" \
  "  The original {ppmoquant.c} was Copyright ©  1989, 1991 by Jef Poskanzer.\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>
#include <argparser.h>
#include <bool.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <uint16_image_RGB_hist.h>
#include <uint16_image_RGB_medcut.h>

#include <ppmoquant.h>

int main(int argc, char* argv[])
  {
    options_t *o = parse_options(argc, argv);

    /* Read the black-bg and white-bg images, ensure that they are compatible: */
    uint16_image_t *bimg = uint16_image_read_pnm_named((o->bimgname == NULL ? "-" : o->bimgname), o->verbose);
    uint16_image_t *wimg = uint16_image_read_pnm_named((o->wimgname == NULL ? "-" : o->wimgname), o->verbose);
    demand(bimg->cols == wimg->cols, "the input images must have the same width");
    demand(bimg->rows == wimg->rows, "the input images must have the same height");
    demand(bimg->maxval == wimg->maxval, "the input images must have the same {maxval}");
      
    /* Read the color map: */
    uint16_image_t *cmap = uint16_image_read_pnm_named(o->cmapname, FALSE);
    int mapmaxval = cmap->maxval;
    int ncolors;
    uint16_image_RGB_hist_vector cm = collect_colors(cmap, &ncolors);
    if (o->verbose) { fprintf(stderr, "using %d colors\n", ncolors); }
    
    /* Make sure that the transparent color is disjoint from the colormap: */
    ppm_pixel_t *trp = find_color(o->transpcolor, cm, ncolors);
    demand(trp == NULL, "\"-transparent\" color should not occur in colormap");
    
    /* Quantize the images: */
    if (o->floyd)
      { uint16_image_quantize_floyds(bimg, wimg, cm, ncolors, o->transpcolor, mapmaxval, o->match); }
    else                                                   
      { uint16_image_quantize_plains(bimg, wimg, cm, ncolors, o->transpcolor, mapmaxval, o->match); }

    uint16_image_write_pnm_named("-", bimg, 0, FALSE);
    exit(0);
  }

options_t *parse_options(int argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
     
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
    
    /* !!! Add the option RGB/YUV and use it. !!! */
    o->match = &rgb_match_range;

    /* Keyword arguments: */
    
    if ( argparser_keyword_present(pp, "-fs") || argparser_keyword_present(pp, "-floyd"))
      { o->floyd = 1; }
    else if (argparser_keyword_present(pp, "-nofs") || argparser_keyword_present(pp, "-nofloyd"))
      { o->floyd = 0; }
      
    argparser_get_keyword(pp, "-map");
    o->cmapname = argparser_get_next(pp);
      
    argparser_get_keyword(pp, "-transparent");
    o->transpcolor = my_parsecolor(pp);

    /* Go to positional arguments: */
    argparser_skip_parsed(pp);

    /* Black-background and white-background images are mandatory: */
    o->bimgname = argparser_get_next(pp);
    o->wimgname = argparser_get_next(pp);
    
    o->verbose = argparser_keyword_present(pp, "-verbose");

    /* Check for spurious args: */
    argparser_finish(pp);
    
    return o;
  }

ppm_pixel_t my_parsecolor(argparser_t *pp)
  { ppm_pixel_t p;
    int c; 
    for (c = 0; c < 3; c++) 
      { p.c[c] = argparser_get_next_int(pp, 0, PNM_FILE_MAX_MAXVAL); }
    return p;
  }

uint16_image_RGB_hist_vector collect_colors(uint16_image_t *map, int *ncolorsp)
  { if (map->cols == 0 || map->rows == 0) { pnm_error("null colormap??"); }
    uint16_image_RGB_hist_vector cm = uint16_image_RGB_hist_build(map->smp, map->chns, map->cols, map->rows, MAXCOLORS, ncolorsp);
    if (cm == NULL) { pnm_error("too many colors in colormap!"); }
    pnm_message("%d colors found in colormap", (*ncolorsp));
    return(cm);
  }

ppm_pixel_t* find_color
  ( ppm_pixel_t p,
    uint16_image_RGB_hist_vector cm, 
    int ncolors
  )
  { register int i;
    register int pr = p.c[0];
    register int pg = p.c[1];
    register int pb = p.c[2];
    for (i = 0; i < ncolors; ++i)
      { ppm_pixel_t *qq = &(cm[i].color);
        int r = (*qq).c[0];
        int g = (*qq).c[1];
        int b = (*qq).c[2];
        if ((r == pr) && (g == pg) && (b == pb)) return (qq);
      }
    return(NULL);
  }

void uint16_image_quantize_floyds
  ( uint16_image_t *bimg, 
    uint16_image_t *wimg, 
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    ppm_pixel_t transpcolor,
    int mapmaxval,
    matchfn_t *match
  )
  {
#define FS_SCALE 1024
#define FS_WTA 7
#define FS_WTB 3
#define FS_WTC 5
#define FS_WTD 1
    fs_errors be, we;
    int chns = bimg->chns; assert(wimg->chns == chns); assert(chns == 3);
    int cols = bimg->cols; assert(wimg->cols == cols);
    int rows = bimg->rows; assert(wimg->rows == rows);
    register uint16_t *bp, *wp;
    int maxval = bimg->maxval;
    long maxcor = 2*maxval;
    long mincor = -1*maxval;
    /* int addtohash = 1; */
    /* colorhash_table cht = ppm_alloccolorhash(); */
    register long br, bg, bb, wr, wg, wb;
    double coef = FS_SCALE*((double)maxval)/((double)mapmaxval);

    floyd_alloc_error_vectors(&be, cols);
    floyd_alloc_error_vectors(&we, wimg->cols);

    /* Initialize Floyd-Steinberg error vectors with randoms in [-1..+1]. */
    srand((int) 46157);
    int col;
    for (col = 0; col < cols + 2; ++col)
      {
        be.thisr[col] = (rand() % (FS_SCALE * 2)) - FS_SCALE;
        be.thisg[col] = (rand() % (FS_SCALE * 2)) - FS_SCALE;
        be.thisb[col] = (rand() % (FS_SCALE * 2)) - FS_SCALE;

        we.thisr[col] = be.thisr[col];
        we.thisg[col] = be.thisg[col];
        we.thisb[col] = be.thisb[col];

      }
      
    int row;
    int fs_direction = 1; /* Alternates between 0 and 1 at each row. */
    for (row = 0; row < rows; ++row)
      { for (col = 0; col < cols + 2; ++col)
          { be.nextr[col] = be.nextg[col] = be.nextb[col] = 0;
            we.nextr[col] = we.nextg[col] = we.nextb[col] = 0;
          }
        int limitcol;
        if (fs_direction)
          { col = 0; limitcol = cols;
            bp = &(bimg->smp[row][0]);
            wp = &(wimg->smp[row][0]);
          }
        else
          { col = cols - 1; limitcol = -1;
            bp = &(bimg->smp[row][col*chns]);
            wp = &(wimg->smp[row][col*chns]);
          }
        do
          { /* Use Floyd-Steinberg errors to adjust actual color. */
            br = bp[0] + be.thisr[col + 1] / FS_SCALE;
            bg = bp[1] + be.thisg[col + 1] / FS_SCALE;
            bb = bp[2] + be.thisb[col + 1] / FS_SCALE;
            
            if (br < mincor) { br = mincor; } if (br > maxcor) { br = maxcor; }
            if (bg < mincor) { bg = mincor; } if (bg > maxcor) { bg = maxcor; }
            if (bb < mincor) { bb = mincor; } if (bb > maxcor) { bb = maxcor; }

            wr = wp[0] + we.thisr[col + 1] / FS_SCALE;
            wg = wp[1] + we.thisg[col + 1] / FS_SCALE;
            wb = wp[2] + we.thisb[col + 1] / FS_SCALE;

            if (wr < mincor) { wr = mincor; } if (wr > maxcor) { wr = maxcor; }
            if (wg < mincor) { wg = mincor; } if (wg > maxcor) { wg = maxcor; }
            if (wb < mincor) { wb = mincor; } if (wb > maxcor) { wb = maxcor; }

            /* Choose replacement color */
            ppm_pixel_t *qq = choose_color
              ( br, bg, bb, wr, wg, wb, maxval,
                cm, ncolors, mapmaxval,
                match
              );
            /* Replace ppm_pixel_t in black image, and propagate error: */
            if (qq == NULL)
              { /* Best match is transparent: */
                bp[0] = transpcolor.c[0];
                bp[1] = transpcolor.c[1];
                bp[2] = transpcolor.c[2];
                floyd_propagate_error_transp
                  ( br, bg, bb,  wr, wg, wb,  
                    maxval,
                    &be, &we, col, fs_direction
                  );
              }
            else
              { /* Best match is opaque: */
                bp[0] = qq->c[0];
                bp[1] = qq->c[1];
                bp[2] = qq->c[2];
                floyd_propagate_error_opaque
                  ( br, bg, bb,  wr, wg, wb,
                    qq, coef,
                    &be, &we, col, fs_direction
                  );
              }

            if (fs_direction)
              { ++col; bp += chns; wp += chns; }
            else
              { --col; bp -= chns; wp -= chns; }
            }
        while (col != limitcol);
        floyd_swap_error_vectors(&be);
        floyd_swap_error_vectors(&we);
        
        fs_direction = ! fs_direction;
      }
    bimg->maxval = mapmaxval;
  }

void floyd_alloc_error_vectors(fs_errors *ep, int cols)
  {
    ep->thisr = (long*)notnull(malloc((cols + 2)*sizeof(long)), "no mem");
    ep->nextr = (long*)notnull(malloc((cols + 2)*sizeof(long)), "no mem");
    ep->thisg = (long*)notnull(malloc((cols + 2)*sizeof(long)), "no mem");
    ep->nextg = (long*)notnull(malloc((cols + 2)*sizeof(long)), "no mem");
    ep->thisb = (long*)notnull(malloc((cols + 2)*sizeof(long)), "no mem");
    ep->nextb = (long*)notnull(malloc((cols + 2)*sizeof(long)), "no mem");
  }

void floyd_propagate_error_opaque
  ( long br, long bg, long bb,
    long wr, long wg, long wb,
    ppm_pixel_t *newp, double coef,
    fs_errors *bep, fs_errors *wep,
    int col,
    int fs_direction
  )
  {
    long err, new;
    if (fs_direction)
      {
        new = (long)(coef*newp->c[0]+0.5);
        err = (br * FS_SCALE - new);
        bep->thisr[col + 2] += (err * FS_WTA) / 16;
        bep->nextr[col    ] += (err * FS_WTB) / 16;
        bep->nextr[col + 1] += (err * FS_WTC) / 16;
        bep->nextr[col + 2] += (err * FS_WTD) / 16;
        err = (wr * FS_SCALE - new);
        wep->thisr[col + 2] += (err * FS_WTA) / 16;
        wep->nextr[col    ] += (err * FS_WTB) / 16;
        wep->nextr[col + 1] += (err * FS_WTC) / 16;
        wep->nextr[col + 2] += (err * FS_WTD) / 16;

        new = (long)(coef*newp->c[1]+0.5);
        err = (bg * FS_SCALE - new);
        bep->thisg[col + 2] += (err * FS_WTA) / 16;
        bep->nextg[col    ] += (err * FS_WTB) / 16;
        bep->nextg[col + 1] += (err * FS_WTC) / 16;
        bep->nextg[col + 2] += (err * FS_WTD) / 16;
        err = (wg * FS_SCALE - new);
        wep->thisg[col + 2] += (err * FS_WTA) / 16;
        wep->nextg[col    ] += (err * FS_WTB) / 16;
        wep->nextg[col + 1] += (err * FS_WTC) / 16;
        wep->nextg[col + 2] += (err * FS_WTD) / 16;

        new = (long)(coef*newp->c[2]+0.5);
        err = (bb * FS_SCALE - new);
        bep->thisb[col + 2] += (err * FS_WTA) / 16;
        bep->nextb[col    ] += (err * FS_WTB) / 16;
        bep->nextb[col + 1] += (err * FS_WTC) / 16;
        bep->nextb[col + 2] += (err * FS_WTD) / 16;
        err = (wb * FS_SCALE - new);
        wep->thisb[col + 2] += (err * FS_WTA) / 16;
        wep->nextb[col    ] += (err * FS_WTB) / 16;
        wep->nextb[col + 1] += (err * FS_WTC) / 16;
        wep->nextb[col + 2] += (err * FS_WTD) / 16;
      }
    else
      {
        new = (long)(coef*newp->c[0]+0.5);
        err = (br * FS_SCALE - new);
        bep->thisr[col    ] += (err * FS_WTA) / 16;
        bep->nextr[col + 2] += (err * FS_WTB) / 16;
        bep->nextr[col + 1] += (err * FS_WTC) / 16;
        bep->nextr[col    ] += (err * FS_WTD) / 16;
        err = (wr * FS_SCALE - new);
        wep->thisr[col    ] += (err * FS_WTA) / 16;
        wep->nextr[col + 2] += (err * FS_WTB) / 16;
        wep->nextr[col + 1] += (err * FS_WTC) / 16;
        wep->nextr[col    ] += (err * FS_WTD) / 16;
        
        new = (long)(coef*newp->c[1]+0.5);
        err = (bg * FS_SCALE - new);
        bep->thisg[col    ] += (err * FS_WTA) / 16;
        bep->nextg[col + 2] += (err * FS_WTB) / 16;
        bep->nextg[col + 1] += (err * FS_WTC) / 16;
        bep->nextg[col    ] += (err * FS_WTD) / 16;
        err = (wg * FS_SCALE - new);
        wep->thisg[col    ] += (err * FS_WTA) / 16;
        wep->nextg[col + 2] += (err * FS_WTB) / 16;
        wep->nextg[col + 1] += (err * FS_WTC) / 16;
        wep->nextg[col    ] += (err * FS_WTD) / 16;
        
        new = (long)(coef*newp->c[2]+0.5);
        err = (bb * FS_SCALE - new);
        bep->thisb[col    ] += (err * FS_WTA) / 16;
        bep->nextb[col + 2] += (err * FS_WTB) / 16;
        bep->nextb[col + 1] += (err * FS_WTC) / 16;
        bep->nextb[col    ] += (err * FS_WTD) / 16;
        err = (wb * FS_SCALE - new);
        wep->thisb[col    ] += (err * FS_WTA) / 16;
        wep->nextb[col + 2] += (err * FS_WTB) / 16;
        wep->nextb[col + 1] += (err * FS_WTC) / 16;
        wep->nextb[col    ] += (err * FS_WTD) / 16;
      }
  }

void  floyd_propagate_error_transp
  ( long br, long bg, long bb,
    long wr, long wg, long wb,
    int maxval,
    fs_errors *bep, fs_errors *wep,
    int col,
    int fs_direction
  )
  {
    long err;
    if (fs_direction)
      {
        err = br * FS_SCALE;
        bep->thisr[col + 2] += (err * FS_WTA) / 16;
        bep->nextr[col    ] += (err * FS_WTB) / 16;
        bep->nextr[col + 1] += (err * FS_WTC) / 16;
        bep->nextr[col + 2] += (err * FS_WTD) / 16;
        err = (wr - maxval) * FS_SCALE;
        wep->thisr[col + 2] += (err * FS_WTA) / 16;
        wep->nextr[col    ] += (err * FS_WTB) / 16;
        wep->nextr[col + 1] += (err * FS_WTC) / 16;
        wep->nextr[col + 2] += (err * FS_WTD) / 16;

        err = bg * FS_SCALE;
        bep->thisg[col + 2] += (err * FS_WTA) / 16;
        bep->nextg[col    ] += (err * FS_WTB) / 16;
        bep->nextg[col + 1] += (err * FS_WTC) / 16;
        bep->nextg[col + 2] += (err * FS_WTD) / 16;
        err = (wg - maxval) * FS_SCALE;
        wep->thisg[col + 2] += (err * FS_WTA) / 16;
        wep->nextg[col    ] += (err * FS_WTB) / 16;
        wep->nextg[col + 1] += (err * FS_WTC) / 16;
        wep->nextg[col + 2] += (err * FS_WTD) / 16;

        err = bb * FS_SCALE;
        bep->thisb[col + 2] += (err * FS_WTA) / 16;
        bep->nextb[col    ] += (err * FS_WTB) / 16;
        bep->nextb[col + 1] += (err * FS_WTC) / 16;
        bep->nextb[col + 2] += (err * FS_WTD) / 16;
        err = (wb - maxval) * FS_SCALE;
        wep->thisb[col + 2] += (err * FS_WTA) / 16;
        wep->nextb[col    ] += (err * FS_WTB) / 16;
        wep->nextb[col + 1] += (err * FS_WTC) / 16;
        wep->nextb[col + 2] += (err * FS_WTD) / 16;
      }
    else
      {
        err = br * FS_SCALE;
        bep->thisr[col    ] += (err * FS_WTA) / 16;
        bep->nextr[col + 2] += (err * FS_WTB) / 16;
        bep->nextr[col + 1] += (err * FS_WTC) / 16;
        bep->nextr[col    ] += (err * FS_WTD) / 16;
        err = (wr - maxval) * FS_SCALE;
        wep->thisr[col    ] += (err * FS_WTA) / 16;
        wep->nextr[col + 2] += (err * FS_WTB) / 16;
        wep->nextr[col + 1] += (err * FS_WTC) / 16;
        wep->nextr[col    ] += (err * FS_WTD) / 16;
        
        err = bg * FS_SCALE;
        bep->thisg[col    ] += (err * FS_WTA) / 16;
        bep->nextg[col + 2] += (err * FS_WTB) / 16;
        bep->nextg[col + 1] += (err * FS_WTC) / 16;
        bep->nextg[col    ] += (err * FS_WTD) / 16;
        err = (wg - maxval) * FS_SCALE;
        wep->thisg[col    ] += (err * FS_WTA) / 16;
        wep->nextg[col + 2] += (err * FS_WTB) / 16;
        wep->nextg[col + 1] += (err * FS_WTC) / 16;
        wep->nextg[col    ] += (err * FS_WTD) / 16;
        
        err = bb * FS_SCALE;
        bep->thisb[col    ] += (err * FS_WTA) / 16;
        bep->nextb[col + 2] += (err * FS_WTB) / 16;
        bep->nextb[col + 1] += (err * FS_WTC) / 16;
        bep->nextb[col    ] += (err * FS_WTD) / 16;
        err = (wb - maxval) * FS_SCALE;
        wep->thisb[col    ] += (err * FS_WTA) / 16;
        wep->nextb[col + 2] += (err * FS_WTB) / 16;
        wep->nextb[col + 1] += (err * FS_WTC) / 16;
        wep->nextb[col    ] += (err * FS_WTD) / 16;
      }
  }
  
void floyd_swap_error_vectors(fs_errors *ep)
  {
    long* t;
    t = ep->thisr; ep->thisr = ep->nextr; ep->nextr = t;
    t = ep->thisg; ep->thisg = ep->nextg; ep->nextg = t;
    t = ep->thisb; ep->thisb = ep->nextb; ep->nextb = t;
  }

void uint16_image_quantize_plains
  ( uint16_image_t *bimg, 
    uint16_image_t *wimg, 
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    ppm_pixel_t transpcolor,
    int mapmaxval,
    matchfn_t *match
  )
  {
    int chns = bimg->chns; assert(wimg->chns == chns); assert(chns == 3);
    int cols = bimg->cols; assert(wimg->cols == cols);
    int rows = bimg->rows; assert(wimg->rows == rows);
    uint16_t *bp, *wp;
    long br, bg, bb, wr, wg, wb;
    int maxval = bimg->maxval;
    
    int row;
    for (row = 0; row < bimg->rows; ++row)
      {
        int col = 0;
        int limitcol = cols;
        bp = bimg->smp[row];
        wp = wimg->smp[row];

        do
          {
            br = bp[0];
            bg = bp[1];
            bb = bp[2];

            wr = wp[0];
            wg = wp[1];
            wb = wp[2];
                
            /* Choose replacement color* */
            ppm_pixel_t *qq = choose_color
              ( br, bg, bb, wr, wg, wb, maxval,
                cm, ncolors, mapmaxval,
                match
              );
            /* Replace ppm_pixel_t in black image, and propagate error: */
            if (qq == NULL)
              { /* Best match is transparent: */
                bp[0] = transpcolor.c[0];
                bp[1] = transpcolor.c[1];
                bp[2] = transpcolor.c[2];
              }
            else
              { /* Best match is opaque: */
                bp[0] = qq->c[0];
                bp[1] = qq->c[1];
                bp[2] = qq->c[2];
              }
            ++col; bp += bimg->chns; wp += wimg->chns;
          }
        while (col != limitcol);
      }
    bimg->maxval = mapmaxval;
  }

ppm_pixel_t* choose_color 
  ( long br, long bg, long bb,
    long wr, long wg, long wb,
    int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval,
    matchfn_t *match
  )
  {
    int ind = match
  ( br, bg, bb, wr, wg, wb, maxval,
      cm, ncolors, mapmaxval
    );
    if (ind == -1) return(NULL);
    else return(&(cm[ind].color));
  }

int rgb_match_range 
  ( long br, long bg, long bb,
    long wr, long wg, long wb,
    int maxval,
    uint16_image_RGB_hist_vector cm, 
    int ncolors,
    int mapmaxval
  )
  { 
#define RGBMATCH_SCALE 4
    register int i;
    register int ind = -1;
    register long tr, tg, tb;
    register long d, dist, newdist;
        
    /* Scale the ppm_pixel_t coordinates to match the map's maxval: */
    double s = RGBMATCH_SCALE * ((double)mapmaxval)/((double)maxval);
    
    long fbr = s*br + 0.5;
    long fbg = s*bg + 0.5;
    long fbb = s*bb + 0.5;

    long fwr = s*wr + 0.5;
    long fwg = s*wg + 0.5;
    long fwb = s*wb + 0.5;
    
    long fmx = s*maxval + 0.5;
    
    /* 
      Initialize dist with the distance between the transparent
      color range and the given color range:
    */
    dist = 0;
    d = fbr; if (d < 0) { d = -d; } if (d > dist) { dist = d; }
    d = fbg; if (d < 0) { d = -d; } if (d > dist) { dist = d; }
    d = fbb; if (d < 0) { d = -d; } if (d > dist) { dist = d; }
    
    d = fwr - fmx; if (d < 0) { d = -d; } if (d > dist) { dist = d; }
    d = fwg - fmx; if (d < 0) { d = -d; } if (d > dist) { dist = d; }
    d = fwb - fmx; if (d < 0) { d = -d; } if (d > dist) { dist = d; }

    /* Now see if an opaque color would work better: */
    for (i = 0; ((dist > 0) && (i < ncolors)); ++i)
      { ppm_pixel_t *qq = &(cm[i].color);
        tr = RGBMATCH_SCALE * (*qq).c[0];
        tg = RGBMATCH_SCALE * (*qq).c[1];
        tb = RGBMATCH_SCALE * (*qq).c[2];
        /* 
          Find the maximum distance newdist between this color and the corners 
          of the given interval. If it is larger than dist, skip it.
        */
        newdist = 0;
        d = tr - fbr; if (d < 0) { d = -d; }
        if (d > newdist) { if (d >= dist) goto skipit; newdist = d; }
        d = tg - fbg; if (d < 0) { d = -d; }
        if (d > newdist) { if (d >= dist) goto skipit; newdist = d; }
        d = tb - fbb; if (d < 0) { d = -d; }
        if (d > newdist) { if (d >= dist) goto skipit; newdist = d; }
            
        d = tr - fwr; if (d < 0) { d = -d; }
        if (d > newdist) { if (d >= dist) goto skipit; newdist = d; }
        d = tg - fwg; if (d < 0) { d = -d; }
        if (d > newdist) { if (d >= dist) goto skipit; newdist = d; }
        d = tb - fwb; if (d < 0) { d = -d; }
        if (d > newdist) { if (d >= dist) goto skipit; newdist = d; }
            
        /* If we got here, newdist is smaller than dist: */
        ind = i; dist = newdist;
      skipit:
        /* OK */;
      }
    return (ind);
  }
