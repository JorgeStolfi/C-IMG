#define PROG_NAME "pnmtopsx"
#define PROG_DESC "convert a PBM/PGM/PPM file to PostScript w/ true scaling"
#define PROG_VERS "2.0"

/* Copyright © 1989 by Jef Poskanzer, modified by J. Stolfi. */
/* See end of file for full copyright and (no)warranty notice. */
/* Last edited on 2023-03-07 18:06:01 by stolfi */

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "  [-scale {SCALE}] [-turn|-noturn] [-center|-nocenter] \\\n" \
  "  [-rle|-runlength] [-dpi {DPI}] \\\n" \
  "  [-width {WIDTH}] [-height {HEIGHT}] \\\n" \
  "  [<] {PNMFILE} \\\n" \
  "  > {PSFILE}"

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  "  Same as pnmtops(1) but obeys the \"-scale\" parameter literally." \
  "OPTIONS\n" \
  "  Same as pnmtops(1), 1994 NetPBM version.\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  pnm(5), psidtopgm(1), pnmtops(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 1989-1991 by Jef Poskanzer.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  mar/2023: int->int32_t, other code cleanup. Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "  nov/2006: rewritten to use leaner PBM libs," \
  " argparser, etc. Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "  oct/1999: Scaling logic modified by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "  mar/1994: Netpbm-1mar1994 release of \"pnmtops\"\n" \
  "\n" \
  "  nov/1993: Option \"-nocenter\" by Wolfgang Stuerzlinger, wrzl@gup.uni-linz.ac.at.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 1989, 1991 by Jef Poskanzer.\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#include <bool.h>
#include <jsfile.h>
#include <jspnm.h>
#include <uint16_image.h>
#include <argparser.h>

/* COMMAND-LINE OPTIONS */

typedef struct options_t
  { bool_t turn_do;
    bool_t turn_ok;
    bool_t rle;
    bool_t center;
    double scale;
    int32_t dpi;
    double width;
    double height;
    char *iname; /* Input file name. */
  } options_t;

typedef struct rgb_pixel_t { uint16_t c[3]; } rgb_pixel_t; 

/* INTERNAL PROTOTYPES */

options_t *parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and packs them as an {options_t}. */

int32_t main(int32_t argc,char** argv);
void put_init(options_t *o, int32_t cols, int32_t rows, int32_t chns, int32_t bps, int32_t pad, char *pname);
void put_item (void);
void put_sample (uint16_t xv);
void put_rest (void);
void rle_put_buffer (void);
void rle_put_item (void);
void rle_put_sample (uint16_t xv);
void rle_flush (void);
void rle_put_rest (void);
char *make_pname_from_iname(char *iname);

/* Was 0.95 -- stolfi */
#define MARGIN 1.0

int32_t main(int32_t argc, char* argv[])
  {
    options_t *o = parse_options(argc, argv);
    
    FILE* ifp = open_read(o->iname, FALSE);
    int32_t rows, cols, chns;
    uint16_t maxval_img;
    bool_t raw,bits;
    pnm_format_t format;
    pnm_read_header(ifp, &cols, &rows, &chns, &maxval_img, &raw, &bits, &format);
    
    uint16_t *smprow = uint16_image_alloc_pixel_row(cols, chns);
    
    /* Figure out bps. */
    int32_t bps = 0;
    uint16_t maxval_eps = 0;
    if (maxval_img <= 1)
      { bps = 1; maxval_eps = 1; }
    else if (maxval_img <= 3)
      { bps = 2; maxval_eps = 3; } 
    else if ( maxval_img <= 15)
      { bps = 4; maxval_eps = 15; }
    else if (maxval_img <= 255)
      { bps = 8; maxval_eps = 255; }
    else
      { fprintf(stderr, "!!maxval %d will be reduced to 255", maxval_img);
        bps = 8; maxval_eps = 255;
      }
    
    /* Compute padding to round {cols * bps} up to the nearest multiple of 8. */
    int32_t padright = (((cols * bps + 7) / 8) * 8 - cols * bps) / bps;

    char *pname = make_pname_from_iname(o->iname);
    
    put_init(o, cols, rows, chns, bps, padright, pname);

    int32_t sprw = cols*chns; /* Samples per row. */
    
    /* Compute scale factor for {maxval_img} to {maxval_eps} recoding: */
    double scale = ((double)maxval_eps)/((double)maxval_img);

    uint16_t* xP;
    for (int32_t row = 0; row < rows; ++row )
      { pnm_read_pixels(ifp, smprow, cols, chns, maxval_img, raw, bits);
        if (maxval_img != maxval_eps)
          { /* Adjust scale from {maxval_img} to {maxval_eps}: */
            xP = smprow;
            for (int32_t k = 0; k < sprw; k++)
              { int32_t remap = (int32_t)floor(((double)(*xP))*scale + 0.5);
                assert(remap <= maxval_eps);
                (*xP) = (int16_t)remap; 
                xP++;
              }
          }
        /* Output row samples, one channel at a time: */
        for (int32_t kc = 0; kc < chns; kc++)
          { /* Output one every {chns} samples, starting with sample {kc}: */
            xP = smprow + kc;
            for (int32_t col = 0; col < cols; col++)
              { if (o->rle)
                  { rle_put_sample(*xP); }
                else
                  { put_sample(*xP); }
                xP += chns;
              }
            /* Pad the data with zeros: */
            for (int32_t col = 0; col < padright; ++col)
              { if (o->rle)
                  { rle_put_sample(0); }
                else
                  { put_sample(0); }
              }
            if (o->rle) { rle_flush(); }
          }
      }

    if (ifp != stdin) { fclose(ifp); } else { fflush(ifp); };

    if (o->rle)
      { rle_put_rest(); }
    else
      { put_rest(); }
    return 0;
  }

static int32_t bitspersample, item, bitsperitem, bitshift, itemsperline, items;
static int32_t rleitem, rlebitsperitem, rlebitshift;
static int32_t repeat, itembuf[128], count, repeatitem, repeatcount;

void put_init(options_t *o, int32_t cols, int32_t rows, int32_t chns, int32_t bps, int32_t pad, char *pname)
  {
    /* Turn? */
    int32_t icols = cols;
    int32_t irows = rows;
    bool_t turn = o->turn_do || (o->turn_ok && (cols > rows));
    if (turn)
      { cols = irows;
        rows = icols;
      }

    /* Figure out size. */
    /* Removed automatic scale adjustment */
    double scols = o->scale * cols;
    double srows = o->scale * rows;
    double llx = (o->center) ? ( o->width - scols ) / 2 : 0;
    double lly = (o->center) ? ( o->height - srows ) / 2 : 0;

    printf( "%%!PS-Adobe-2.0 EPSF-2.0\n" );
    printf( "%%%%Creator: pnmtopsx\n" );
    printf( "%%%%Title: %s.ps\n", pname );
    printf( "%%%%Pages: 1\n" );
    printf(
        "%%%%BoundingBox: %d %d %d %d\n",
        (int32_t)floor(llx), (int32_t)floor(lly),
        (int32_t)ceil(llx + scols), (int32_t)ceil(lly + srows) );
    printf( "%%%%EndComments\n" );
    if ( o->rle )
        {
        printf( "/rlestr1 1 string def\n" );
        printf( "/readrlestring {\n" );                         /* s -- nr */
        printf( "  /rlestr exch def\n" );                       /* - */
        printf( "  currentfile rlestr1 readhexstring pop\n" );  /* s1 */
        printf( "  0 get\n" );                                  /* c */
        printf( "  dup 127 le {\n" );                           /* c */
        printf( "    currentfile rlestr 0\n" );                 /* c f s 0 */
        printf( "    4 3 roll\n" );                             /* f s 0 c */
        printf( "    1 add  getinterval\n" );                   /* f s */
        printf( "    readhexstring pop\n" );                    /* s */
        printf( "    length\n" );                               /* nr */
        printf( "  } {\n" );                                    /* c */
        printf( "    256 exch sub dup\n" );                     /* n n */
        printf( "    currentfile rlestr1 readhexstring pop\n" );/* n n s1 */
        printf( "    0 get\n" );                                /* n n c */
        printf( "    exch 0 exch 1 exch 1 sub {\n" );           /* n c 0 1 n-1*/
        printf( "      rlestr exch 2 index put\n" );
        printf( "    } for\n" );                                /* n c */
        printf( "    pop\n" );                                  /* nr */
        printf( "  } ifelse\n" );                               /* nr */
        printf( "} bind def\n" );
        printf( "/readstring {\n" );                            /* s -- s */
        printf( "  dup length 0 {\n" );                         /* s l 0 */
        printf( "    3 copy exch\n" );                          /* s l n s n l*/
        printf( "    1 index sub\n" );                          /* s l n s n r*/
        printf( "    getinterval\n" );                          /* s l n ss */
        printf( "    readrlestring\n" );                        /* s l n nr */
        printf( "    add\n" );                                  /* s l n */
        printf( "    2 copy le { exit } if\n" );                /* s l n */
        printf( "  } loop\n" );                                 /* s l l */
        printf( "  pop pop\n" );                                /* s */
        printf( "} bind def\n" );
        }
    else
        {
        printf( "/readstring {\n" );                            /* s -- s */
        printf( "  currentfile exch readhexstring pop\n" );
        printf( "} bind def\n" );
        }

    if (chns == 3)
        {
        printf( "/rpicstr %d string def\n", ( icols + pad ) * bps / 8 );
        printf( "/gpicstr %d string def\n", ( icols + pad ) * bps / 8 );
        printf( "/bpicstr %d string def\n", ( icols + pad ) * bps / 8 );
        }
    else
        printf( "/picstr %d string def\n", ( icols + pad ) * bps / 8 );

    printf( "%%%%EndProlog\n" );
    printf( "%%%%Page: 1 1\n" );
    printf( "gsave\n" );
    printf( "%g %g translate\n", llx, lly );
    printf( "%g %g scale\n", scols, srows );
    if ( turn ) printf( "0.5 0.5 translate  90 rotate  -0.5 -0.5 translate\n" );
    printf( "%d %d %d\n", icols, irows, bps );
    printf( "[ %d 0 0 -%d 0 %d ]\n", icols, irows, irows );

    if (chns == 3)
        {
        printf( "{ rpicstr readstring }\n" );
        printf( "{ gpicstr readstring }\n" );
        printf( "{ bpicstr readstring }\n" );
        printf( "true 3\n" );
        printf( "colorimage\n" );
        pnm_message( "writing color PostScript..." );
        }
    else
        {
        printf( "{ picstr readstring }\n" );
        printf( "image\n" );
        }

    bitspersample = bps;
    itemsperline = items = 0;
    if ( o->rle )
        {
        rleitem = 0;
        rlebitsperitem = 0;
        rlebitshift = 8 - bitspersample;
        repeat = 1;
        count = 0;
        }
    else
        {
        item = 0;
        bitsperitem = 0;
        bitshift = 8 - bitspersample;
        }
  }

void put_item(void)
  {
    char* hexits = "0123456789abcdef";

    if ( itemsperline == 30 )
        {
        putchar( '\n' );
        itemsperline = 0;
        }
    putchar( hexits[item >> 4] );
    putchar( hexits[item & 15] );
    ++itemsperline;
    ++items;
    item = 0;
    bitsperitem = 0;
    bitshift = 8 - bitspersample;
  }

void put_sample( uint16_t xv )
  {
    if ( bitsperitem == 8 )
        put_item();
    item += xv << bitshift;
    bitsperitem += bitspersample;
    bitshift -= bitspersample;
  }

void put_rest()
  {
    if ( bitsperitem > 0 )
        put_item();
    printf( "\n" );
    printf( "grestore\n" );
    printf( "showpage\n" );
    printf( "%%%%Trailer\n" );
  }

void rle_put_buffer()
  {
    if ( repeat )
        {
        item = 256 - count;
        put_item();
        item = repeatitem;
        put_item();
        }
    else
        {
        item = count - 1;
        put_item();
        for (int32_t i = 0; i < count; ++i )
            {
            item = itembuf[i];
            put_item();
            }
        }
    repeat = 1;
    count = 0;
  }

void rle_put_item()
  {
    if ( count == 128 )
        rle_put_buffer();

    if ( repeat && count == 0 )
        { /* Still initializing a repeat buf. */
        itembuf[count] = repeatitem = rleitem;
        ++count;
        }
    else if ( repeat )
        { /* Repeating - watch for end of run. */
        if ( rleitem == repeatitem )
            { /* Run continues. */
            itembuf[count] = rleitem;
            ++count;
            }
        else
            { /* Run ended - is it long enough to dump? */
            if ( count > 2 )
                { /* Yes, dump a repeat-mode buffer and start a new one. */
                rle_put_buffer();
                itembuf[count] = repeatitem = rleitem;
                ++count;
                }
            else
                { /* Not long enough - convert to non-repeat mode. */
                repeat = 0;
                itembuf[count] = repeatitem = rleitem;
                ++count;
                repeatcount = 1;
                }
            }
        }
    else
        { /* Not repeating - watch for a run worth repeating. */
        if ( rleitem == repeatitem )
            { /* Possible run continues. */
            ++repeatcount;
            if ( repeatcount > 3 )
                { /* Long enough - dump non-repeat part and start repeat. */
                count = count - ( repeatcount - 1 );
                rle_put_buffer();
                count = repeatcount;
                for (int32_t  i = 0; i < count; ++i )
                    itembuf[i] = rleitem;
                }
            else
                { /* Not long enough yet - continue as non-repeat buf. */
                itembuf[count] = rleitem;
                ++count;
                }
            }
        else
            { /* Broken run. */
            itembuf[count] = repeatitem = rleitem;
            ++count;
            repeatcount = 1;
            }
        }

    rleitem = 0;
    rlebitsperitem = 0;
    rlebitshift = 8 - bitspersample;
  }

void rle_put_sample(uint16_t xv)
  {
    if ( rlebitsperitem == 8 )
        rle_put_item();
    rleitem += xv << rlebitshift;
    rlebitsperitem += bitspersample;
    rlebitshift -= bitspersample;
  }

void rle_flush()
  {
    if ( rlebitsperitem > 0 )
        rle_put_item();
    if ( count > 0 )
        rle_put_buffer();
  }

void rle_put_rest()
  {
    rle_flush();
    printf( "\n" );
    printf( "grestore\n" );
    printf( "showpage\n" );
    printf( "%%%%Trailer\n" );
  } 


options_t *parse_options(int32_t argc, char **argv)
  {
    /* Initialize argument parser: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    /* Allocate the program's option record: */
    options_t *o = (options_t *)malloc(sizeof(options_t)); 

    if (argparser_keyword_present(pp, "-scale"))
      { o->scale = argparser_get_next_double(pp, 1.0e-6, 1.0e+6); }
    else
      { o->scale = 1.0; }

    if (argparser_keyword_present(pp, "-turn"))
      { o->turn_do = TRUE;
        o->turn_ok = TRUE;
      }
    else if (argparser_keyword_present(pp, "-noturn"))
      { o->turn_do = FALSE;
        o->turn_ok = FALSE;
      }
    else
      { o->turn_do = FALSE;
        o->turn_ok = TRUE;
      }

    if (argparser_keyword_present(pp, "-center"))
      { o->center = TRUE; }
    else if (argparser_keyword_present(pp, "-nocenter"))
      { o->center = FALSE; }
    else
      { o->center = TRUE; }

    o->rle = 
       argparser_keyword_present(pp, "-rle") || 
       argparser_keyword_present(pp, "-runlength");

    if (argparser_keyword_present(pp, "-dpi"))
      { o->dpi = (int32_t)argparser_get_next_int(pp, 1, INT32_MAX); }
    else
      { o->dpi = 300; /* LaserWriter default. */ }

    if (argparser_keyword_present(pp, "-width"))
      { o->width = argparser_get_next_double(pp, 1.0e-6, 1.0e+6) * 72.0; }
    else
      { o->width = 612; /* LaserWriter default. */ }
   
    if (argparser_keyword_present(pp, "-height"))
      { o->height = argparser_get_next_double(pp, 1.0e-6, 1.0e+6) * 72.0; }
    else
      { o->height = 762; /* LaserWriter default. */ }

    /* Get positional alrguments: */
    argparser_skip_parsed(pp);
   
    if (argparser_next(pp) != NULL)
      { o->iname = argparser_get_next(pp);
        if (strlen(o->iname) == 0)
          { argparser_error(pp, "empty input filename"); }
        else if ((strcmp(o->iname, ".") == 0) || (strcmp(o->iname, "..") == 0))
          { argparser_error(pp, "input file is a directory"); }
      }
    else
      { o->iname = "-"; }

    /* Check for spurious arguments: */
    argparser_finish(pp);
    
    return o;
  }

char *make_pname_from_iname(char *iname)
  { char *pname;
    if (strcmp(iname, "-") == 0) 
      { /* Input is {stdin}, set a valid name for output; */
        pname = "noname";
      }
    else
      { /* Get a copy of the file name: */
        pname = (char *)notnull(malloc(strlen(iname)+1), "no mem");
        strcpy(pname, iname);
        /* Remove the filename's extension, if any: */ 
        char *cp = strrchr(pname+1, '.');
        if (cp != NULL) { *cp = '\0'; }
      }
    return pname;
  }
