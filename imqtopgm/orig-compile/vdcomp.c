/* Last edited on 2023-02-07 22:10:29 by stolfi */
/********************************************************************/
/*  Image Decompression Program - C Version for PC, VAX, Unix       */
/*  and Macintosh systems.                                          */
/*                                                                  */
/*  Decompresses images using Kris Becker's subroutine DECOMP.C     */
/*  which is included in this program in a shortened version.       */
/*                                                                  */
/*  Reads a variable length compressed PDS image and outputs a      */
/*  fixed length uncompressed image file in PDS format with         */
/*  labels, image histogram, engineering table, line header table   */
/*  and an image with PGM labels.   */
/*                                                                  */
/********************************************************************/
/*                                                                  */
/*  Use the following commands to compile and link to produce an    */
/*  executable file:                                                */
/*                                                                  */
/*  On an IBM PC (using Microsoft C Version 5.x)                    */
/*                                                                  */
/*    cl /c vdcomp.c                                                */
/*    link  vdcomp/stack:10000;                                     */
/*                                                                  */
/*  On a VAX:                                                       */
/*                                                                  */
/*    cc   vdcomp                                                   */
/*    $define lnk$library sys$library:vaxcrtl.olb                   */
/*    link vdcomp                                                   */
/*                                                                  */
/*  On a Unix host (Sun, Masscomp)                                  */
/*                                                                  */
/*    cc -o vdcomp vdcomp.c                                         */
/*                                                                  */
/*  On a Macintosh (using Lightspeed C)                             */
/*                                                                  */
/*    link with the following libraries:                            */
/*    stdio, storage, strings, unix and MacTraps, with MacTraps     */
/*    and vdcomp in a separate segment.                             */
/*                                                                  */
/********************************************************************/
/*                                                                  */
/*  Use the following command to run the program:                   */
/*                                                                  */
/*    VDCOMP < [infile] > [outfile]                */
/*                                                                  */
/* HIST                                                             */
/*  feb2023 Greatly simplified the program to write PGM only. */
/*  DEC89 Modified program to handle both Voyager and Viking images.*/
/*  OCT89 Converted Voyager decompression program to handle Viking  */
/*  compressed images.  Changed obuf to 'unsigned' to simplify      */
/*  computation of checksum.                                        */
/*  AUG89 Added code to get command line arguments for filenames    */
/*  and output format; routines to free memory used by the Huffman  */
/*  tree); fixed the SFDU label output length; and modified the     */
/*  I/O routines so that the open for Host type 2 uses binary I/O.  */
/*  JUN89 Fixed READVAR, to get length on 16-bit unswapped hosts.   */
/*  JUL88 C driver to decompress standard Voyager Compressed images */
/*  by Mike Martin 1989/12/02                                       */
/*                                                                  */
/*  Inputs   - Input file to be decompressed.                       */
/*                                                                  */
/*  Outputs  - Output file containing decompressed image.           */
/*                                                                  */
/********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <imq.h>
#include <codetree.h>
                                    /* pc i/o defines               */
#define O_RDONLY         0x0000     /* open for reading only        */
#define O_BINARY         0x8000     /* file mode is binary          */

                                    /* vax i/o defines              */
#define RECORD_TYPE      "rfm=fix"  /* VAX fixed length output      */
#define CTX              "ctx=bin"  /* no translation of \n         */
#define FOP          "fop=cif,sup"  /* file processing ops          */

typedef struct leaf
  {
   struct leaf *right;
   int16_t dn;
   struct leaf *left;
  } NODE;

/*************************************************************************
 Declare the tree pointer. This pointer will hold the root of the tree
 once the tree is created by the accompanying routine huff_tree.
**************************************************************************/

typedef uint8_t byte;

NODE *tree;
int16_t znull = -1;  /* Node {nd} is internal iff {nd->dn} is {znull} */
    
/* subroutine definitions                                           */

void no_labels(void);
uint32_t swap_uint32_t(uint32_t inval);
void decompress(uint32_t line, byte *ibuf, byte *obuf, uint32_t *nin, uint32_t *nout);
void decmpinit(uint32_t *hist);
void free_tree(uint32_t *nfreed);

NODE *huff_tree(uint32_t hist[]);
void dcmprs(uint32_t line, byte *ibuf, byte *obuf, uint32_t *nin, uint32_t *nout, NODE *root);
void sort_freq(uint32_t *freq_list, NODE **node_list, int32_t num_freq);
NODE *new_node(int16_t value);
uint32_t free_node(NODE *pnode, uint32_t total_free);
void print_queue(char *tag, uint32_t *fp, NODE **np, uint32_t num_freq);
void print_queue_entry(uint32_t fp, NODE *np);
void print_node(NODE* np);
void print_codes(NODE* np, char *pref);

/* global variables                                                 */

FILE *infile;
FILE *outfile;
uint32_t                record_bytes, max_lines;
uint32_t                line_samples;
uint32_t               label_checksum = 0L, checksum = 0L;

int main(int argc, char **argv)
{
int32_t buf_size = 2048;
byte ibuf[buf_size],obuf[buf_size];

uint32_t  length,total_bytes,line;
uint32_t  long_length;
uint32_t  count;
int32_t i;

infile = stdin;
outfile = stdout;

/*********************************************************************/
/*                                                                   */
/* read and edit compressed file labels                              */
/*                                                                   */
/*********************************************************************/

   no_labels();

  if (record_bytes == 836)  /* set up values for image sizes */
     {
      max_lines    =  800;
      line_samples =  800;
     }
   else
     {
      max_lines    = 1056;         
      line_samples = 1204;
     }

  fprintf(outfile, "P5\n");
  fprintf(outfile, "%u %u\n", line_samples, max_lines);
  fprintf(outfile, "255\n");
  fflush(outfile);
/*********************************************************************/
/*                                                                   */
/* process the image histogram                                       */
/*                                                                   */
/*********************************************************************/

/* need to know record_bytes,hist_count,hist_item_type,item_count.*/
   total_bytes = 0;
   length = imq_read_var_length_record(infile, buf_size, ibuf);
   total_bytes = total_bytes + length;

   if (record_bytes == 836) /* read one more time for Voyager image */
     {
      length = imq_read_var_length_record(infile, buf_size-record_bytes, ibuf+record_bytes);
      total_bytes = total_bytes + length;
     }

/*********************************************************************/
/*                                                                   */
/* process the encoding histogram                                    */
/* don't have to byte-swap because DECOMP.C does it for us           */
/*                                                                   */
/*********************************************************************/
  uint32_t nhist = 511;
  uint32_t hist[nhist];
  uint32_t hist_bytes = nhist*4;
  
  if (record_bytes == 836)
    {
     length = imq_read_var_length_record(infile, hist_bytes, (byte *)hist);
     length = imq_read_var_length_record(infile, hist_bytes-836, (byte *)hist+836);
     length = imq_read_var_length_record(infile, hist_bytes-1672, (byte *)hist+1672);
    }
  else
    {
     length = imq_read_var_length_record(infile, hist_bytes, (byte *)hist);
     length = imq_read_var_length_record(infile, hist_bytes-1204, (byte *)hist+1204);
    }

/*********************************************************************/
/*                                                                   */
/* process the engineering summary                                   */
/*                                                                   */
/*********************************************************************/

   total_bytes = 0;
   length = imq_read_var_length_record(infile, buf_size, ibuf);

/*********************************************************************/
/*                                                                   */
/* process the line header table                                     */
/*                                                                   */
/*********************************************************************/

  if (record_bytes == 1204)
   {
    long_length = 0L;
    for (i=0;i<1056;i++)
      {
       length = imq_read_var_length_record(infile, buf_size, ibuf);
      }
   }
/*********************************************************************/
/*                                                                   */
/* initialize the decompression                                      */
/*                                                                   */
/*********************************************************************/

    fprintf(stderr, "vdcomp: Initializing decompression routine...\n");
	decmpinit(hist);

/*********************************************************************/
/*                                                                   */
/* decompress the image                                              */
/*                                                                   */
/*********************************************************************/

	fprintf(stderr, "vdcomp: Decompressing data...\n");
    line=0;
    do
      {
       length = imq_read_var_length_record(infile, buf_size, ibuf);
       if (length <= 0) break;
       long_length = (uint32_t)length;
       line += 1;
       decompress(line, ibuf, obuf,&long_length, &record_bytes);
       
      count = (uint16_t)fwrite(obuf, line_samples, 1, outfile);
      if (count != 1)
        { fprintf(stderr, "vdcomp: Error writing output file.  Aborting program.\n");
          fprintf(stderr, "vdcomp: Check disk space or for duplicate file name on VAX.\n");
          exit(1);
        }

       if (record_bytes == 1204) /* do checksum for viking */
         for (i=0;i<record_bytes;i++) checksum += (uint32_t)obuf[i];

       if (line % 100 == 0) fprintf(stderr, "vdcomp:   line %d\n",line);
      } while (length > 0 && line < max_lines);

      if  (record_bytes == 1204) /* print checksum for viking */
      fprintf(stderr, "vdcomp:   image label checksum = %u computed checksum = %u\n",
             label_checksum,checksum);

 fprintf(stderr, "\n");
 free_tree(&long_length);
 fclose(infile);
 fclose(outfile);
}

/*********************************************************************/
/*                                                                   */
/* subroutine no_labels - parse {checksum} and {record_bytes},  */
/* don't write anything  */
/*                                                                   */
/*********************************************************************/

void no_labels(void)
{
int32_t buf_size = 2048;
byte ibuf[buf_size];
uint32_t length;
int32_t i;

do
  {
   length = imq_read_var_length_record(infile, buf_size, ibuf);
   /*****************************************************************/
   /* find the checksum and store in label_checksum                 */
   /*****************************************************************/
   if ((i = strncmp((char*)ibuf," CHECKSUM",9)) == 0)
     { 
       ibuf[length]   = '\0';
       label_checksum = (uint32_t)atol((char*)ibuf+35);
     }

   else if ((i = strncmp((char*)ibuf,"RECORD_BYTES",12)) == 0)
   /*****************************************************************/
   /* get the record_bytes value                                    */
   /*****************************************************************/
     {
      sscanf((char*)ibuf+35,"%d",&record_bytes);
      if (record_bytes != 836) record_bytes = 1204;
     }

   /*****************************************************************/
   /* read to the end of the PDS labels                             */
   /*****************************************************************/
   if ((i = strncmp((char*)ibuf,"END",3)) == 0 && length == 3) break;
  } while (length > 0);

}

uint32_t swap_uint32_t(uint32_t inval)  /* swap 4 byte integer                       */
{
union /* this union is used to swap 16 and 32 bit integers          */
  {
   char  ichar[4];
   uint16_t slen;
   uint32_t  llen;
  } onion;
  char   temp;

  /* byte swap the input field                                      */
  onion.llen   = inval;
  temp   = onion.ichar[0];
  onion.ichar[0]=onion.ichar[3];
  onion.ichar[3]=temp;
  temp   = onion.ichar[1];
  onion.ichar[1]=onion.ichar[2];
  onion.ichar[2]=temp;
  return (onion.llen);
}

void decompress(uint32_t line, byte *ibuf, byte *obuf, uint32_t *nin, uint32_t *nout)
/****************************************************************************
*_TITLE decompress - decompresses image lines stored in compressed format   *
*_ARGS  TYPE       NAME      I/O        DESCRIPTION                         */
   /*      char       *ibuf;  I         Compressed data buffer              */
   /*      char       *obuf;  O         Decompressed image line             */
   /*      int32_t   *nin;    I         Number of bytes on input buffer     */
        /* nout; I         Number of bytes in output buffer    */

  {

/*************************************************************************
  This routine is fairly simple as it's only function is to call the
  routine dcmprs.
**************************************************************************/

    dcmprs(line,ibuf,obuf,nin,nout,tree);

    return;
  }

void decmpinit(uint32_t *hist)
/***************************************************************************
*_TITLE decmpinit - initializes the Huffman tree                           *
*_ARGS  TYPE       NAME      I/O        DESCRIPTION                        */
      /* I     hist;        First-difference histogram.        */

{
  extern NODE *tree;          /* Huffman tree root pointer */

  /* Specify the calling function to initialize the tree */

/****************************************************************************
  Simply call the huff_tree routine and return.
*****************************************************************************/

  tree = huff_tree(hist);

  return;
 }

NODE *huff_tree(uint32_t hist[])
/****************************************************************************
*_TITLE huff_tree - constructs the Huffman tree; returns pointer to root    *
*_ARGS  TYPE          NAME        I/O   DESCRIPTION                         */
      /*     int32_t     *hist;   I    First difference histogram          */

  {
  /*  Local variables used */
    uint32_t freq_list[512];      /* Histogram frequency list */
    NODE **node_list;             /* DN pointer array list */

    register uint32_t *fp;        /* Frequency list pointer */
    register NODE **np;           /* Node list pointer */

    register uint32_t num_freq;   /* Number non-zero frequencies in histogram */

    register uint32_t num_nodes; /* Counter for DN initialization */
    register uint32_t cnt;       /* Miscellaneous counter */

    register NODE *temp;          /* Temporary node pointer */

/***************************************************************************
  Allocate the array of nodes from memory and initialize these with numbers
  corresponding with the frequency list.  There are only 511 possible
  permutations of first difference histograms.  There are 512 allocated
  here to adhere to the FORTRAN version.
****************************************************************************/

   fp = freq_list;
   node_list = (NODE **) malloc(sizeof(temp)*512);
   if (node_list == NULL)
    {
      fprintf(stderr, "** out of memory in huff_tree!\n");
      exit(1);
    }
   np = node_list;

   for (num_nodes=1, cnt=512 ; cnt-- ; num_nodes++)
     {
/**************************************************************************
    The following code has been added to standardize the VAX byte order
    for the "int32_t" type.  This code is intended to make the routine
    as machine independant as possible.
***************************************************************************/
        byte *cp = (byte *) hist++;
        uint32_t j;
        int16_t i;
        for (i=4 ; --i >= 0 ; j = (j << 8) | *(cp+i));

/* Now make the assignment */
        *fp++ = j;
        temp = new_node((int16_t)num_nodes);
        *np++ = temp;
     }

     (*--fp) = 0;         /* Ensure the last element is zeroed out.  */

/***************************************************************************
  Now, sort the frequency list and eliminate all frequencies of zero.
****************************************************************************/

  int32_t debug = 0;

  num_freq = 512;
  sort_freq(freq_list,node_list,num_freq);

  fp = freq_list;
  np = node_list;

  for (num_freq=512 ; (*fp) == 0 && (num_freq) ; fp++, np++, num_freq--);
  if (debug) { print_queue("S", fp,np,num_freq); }

/***************************************************************************
  Now create the tree.  Note that if there is only one difference value,
  it is returned as the root.  On each interation, a new node is created
  and the least frequently occurring difference is assigned to the right
  pointer and the next least frequency to the left pointer.  The node
  assigned to the left pointer now becomes the combination of the two
  nodes and it's frequency is the sum of the two combining nodes.
****************************************************************************/

  for (temp=(*np) ; (num_freq--) > 1 ; )
    {
        temp = new_node(znull);
        temp->right = (*np++);
        temp->left = (*np);
        *np = temp;
        *(fp+1) = *(fp+1) + *fp;
        *fp++ = 0;
        if (debug) { print_queue("M", fp,np,num_freq); }
        sort_freq(fp,np,num_freq);
        if (debug) { print_queue("\nS", fp,np,num_freq); }
        
    }
  if (debug) 
    { fprintf(stderr, "vdcomp: bit codes:\n");
      print_codes(temp, "");
      fprintf(stderr, "\n");
    }
  return temp;
 }
 
void print_codes(NODE* np, char *pref)
  {
    if (np == NULL)
      { fprintf(stderr, "%12s (%s)\n", "NULL", pref); }
    else if (np->dn == znull)
      { int32_t n = (int32_t)strlen(pref);
        char prefx[n+2];
        strcpy(prefx, pref);
        prefx[n+1] = 0;
        prefx[n] = '1';
        print_codes(np->left, prefx);
        prefx[n] = '0';
        print_codes(np->right, prefx);
      }
    else
      { fprintf(stderr, "%+12d (%s)\n", - np->dn + 256, pref); }
  }
  
void print_queue(char *tag, uint32_t *fp, NODE **np, uint32_t num_freq)
  {
    fprintf(stderr, "%s ", tag);
    for(int32_t i = 0; i < num_freq; i++)
      { if ((i < 3) || (i == num_freq-1))
          { if (i > 0) { fprintf(stderr, " "); }
            print_queue_entry(*(fp+i),*(np+i));
          }
        else if ((num_freq > 4) && (i == 3))
          { fprintf(stderr, " ... "); }
      }
    if ((*np)->dn == znull)
      { /* Print children of first node: */
        fprintf(stderr, " ");
        print_node(*np);
        fprintf(stderr, "=(");
        print_node((*np)->left);
        fprintf(stderr, ",");
        print_node((*np)->right);
        fprintf(stderr, ")");
      }
    fprintf(stderr, "\n");
  }

void print_queue_entry(uint32_t fp, NODE *np)
  { fprintf(stderr, "[");
    print_node(np);
    fprintf(stderr, ":%u", fp);
    fprintf(stderr, "]");
  }
    
void print_node(NODE* np)
  { 
    if (np->dn == znull)
      { fprintf(stderr, "@%04lu", ((uint64_t)np) % 10000); }
    else
      { fprintf(stderr, "%+d", np->dn - 256); }
  }

NODE *new_node(int16_t value)
/****************************************************************************
*_TITLE new_node - allocates a NODE structure and returns a pointer to it   *
*_ARGS  TYPE        NAME        I/O     DESCRIPTION                         */
/*                  value;    I      Value to assign to DN field         */

  {
    NODE *temp;         /* Pointer to the memory block */

/***************************************************************************
  Allocate the memory and intialize the fields.
****************************************************************************/

  temp = (NODE *) malloc(sizeof(NODE));

  if (temp != NULL)
    {
      temp->right = NULL;
      temp->dn = value;
      temp->left = NULL;
    }
  else
    {
       fprintf(stderr, "** out of memory in new_node!\n");
       exit(1);
    }

   return temp;
  }

void sort_freq(uint32_t *freq_list, NODE **node_list, int32_t num_freq)
/****************************************************************************
*_TITLE sort_freq - sorts frequency and node lists in increasing freq. order*
*_ARGS  TYPE       NAME            I/O  DESCRIPTION                         */
       /*   freq_list;   I   Pointer to frequency list           */
       /*   node_list;   I   Pointer to array of node pointers   */
       /*   num_freq;    I   Number of values in freq list       */

  {
    /* Local Variables */
    register uint32_t *i;       /* primary pointer into freq_list */
    register uint32_t *j;       /* secondary pointer into freq_list */

    register NODE **k;          /* primary pointer to node_list */
    register NODE **l;          /* secondary pointer into node_list */

    int32_t temp1;             /* temporary storage for freq_list */
    NODE *temp2;                /* temporary storage for node_list */

    register int32_t cnt;      /* count of list elements */

/************************************************************************
  Save the current element - starting with the second - in temporary
  storage.  Compare with all elements in first part of list moving
  each up one element until the element is larger.  Insert current
  element at this point in list.
*************************************************************************/

   if (num_freq <= 0) return;      /* If no elements or invalid, return */

   for (i=freq_list, k=node_list, cnt=num_freq ; --cnt ; *j=temp1, *l=temp2)
     {
        temp1 = *(++i);
        temp2 = *(++k);

        for (j = i, l = k ;  *(j-1) > temp1 ; )
          {
            *j = *(j-1);
            *l = *(l-1);
            j--;
            l--;
            if ( j <= freq_list) break;
          }

     }
  return;
  }

void dcmprs(uint32_t line, byte *ibuf, byte *obuf, uint32_t *nin, uint32_t *nout, NODE *root)
/****************************************************************************
*_TITLE dcmprs - decompresses Huffman coded compressed image lines          *
*_ARGS  TYPE       NAME       I/O       DESCRIPTION                         */
     /*    char       *ibuf;   I        Compressed data buffer              */
     /*    char       *obuf;   O        Decompressed image line             */
     /*    int32_t   *nin;     I        Number of bytes on input buffer     */
     /*    int32_t   *nout;    I        Number of bytes in output buffer    */
     /*    NODE       *root;   I        Huffman coded tree                  */

  {
    /* Local Variables */
    register NODE *ptr = root;        /* pointer to position in tree */
    register byte test;      /* test byte for bit set */
    register byte idn;       /* input compressed byte */
    register byte odn;                /* last pixel value decompressed */
    register byte ndn;                /* next pixel value */

    byte *ilim = ibuf + *nin;         /* end of compressed bytes */
    byte *olim = obuf + *nout;        /* end of output buffer */
    
    int32_t debug = 0;
    int32_t debug2 = 0;

/**************************************************************************
  Check for valid input values for nin, nout and make initial assignments.
***************************************************************************/

    if (ilim > ibuf && olim > obuf)
       odn = *obuf++ = *ibuf++;
    else
       {
           fprintf(stderr, "** invalid byte count in dcmprs!\n");
           exit(1);
       }

/**************************************************************************
  Decompress the input buffer.  Assign the first byte to the working
  variable, idn.  An arithmatic and (&) is performed using the variable
  'test' that is bit shifted to the right.  If the result is 0, then
  go to right else go to left.
***************************************************************************/

    uint32_t nin1 = (*nin)-1;
    if (debug && (line == 1)) 
      { fprintf(stderr, "vdcomp: decompressing line %u", line); 
        fprintf(stderr, "  read %u bytes\n", *nin);
        fprintf(stderr, "  expecting %u pixels...\n", *nout);
        fprintf(stderr, "first few Huffman code bits = ");
        codetree_print_bits(stderr, (nin1 < 3 ? nin1 : 3), ibuf, " ");
        if (nin1 > 3) { fprintf(stderr, " ..."); }
        fprintf(stderr, "\n");
        
      }
    for (idn=(*ibuf) ; ibuf < ilim  ; idn =(*++ibuf))
     {
        for (test=0x80 ; test ; test >>= 1)
           {
            ptr = (test & idn) ? ptr->left : ptr->right;
            if (debug2) 
              { fprintf(stderr, "   %02x %02x %c\n", test, idn, ((test & idn) ? '1' : '0')); }

            if (ptr->dn != -1)
              {
                if (obuf >= olim) return;
                assert((ptr->dn >= 1) && (ptr->dn <= 511));
                int32_t di = ptr->dn - 256;
                ndn = (byte)((int32_t)odn - di);
                if (debug && (line == 1)) 
                  { fprintf(stderr, " %4u %+4d = %4u\n", odn, -di, ndn); }
                *obuf++ = ndn;
                odn = ndn;
                ptr = root;
              }
          }
     }
   if (debug && (line == 1)) { fprintf(stderr, "\n"); }
   return;
  }

void free_tree(uint32_t *nfreed)
/****************************************************************************
*_TITLE free_tree - free memory of all allocated nodes                      *
*_ARGS  TYPE       NAME       I/O        DESCRIPTION                        */
      /*   nfreed;  O        Return of total count of nodes     *
*                                        freed.                             */

/*
*_DESCR This routine is supplied to the programmer to free up all the       *
*       allocated memory required to build the huffman tree.  The count     *
*       of the nodes freed is returned in the parameter 'nfreed'.  The      *
*       purpose of the routine is so if the user wishes to decompress more  *
*       than one file per run, the program will not keep allocating new     *
*       memory without first deallocating all previous nodes associated     *
*       with the previous file decompression.                               *

*_HIST  16-AUG-89 Kris Becker   USGS, Flagstaff Original Version            *
*_END                                                                       *
****************************************************************************/

{
	int32_t total_free = 0;

/****************************************************************************
  Simply call the free_node routine and return the result.
*****************************************************************************/

	*nfreed = free_node(tree,total_free);

	return;
}

uint32_t free_node(NODE *pnode, uint32_t total_free)
/***************************************************************************
*_TITLE free_node - deallocates an allocated NODE pointer
*_ARGS  TYPE     NAME          I/O   DESCRIPTION                           */
        /* NODE     *pnode;       I  Pointer to node to free               */
      /*   int32_t total_free;   I  Total number of freed nodes           */

/*
*_DESCR  free_node will check both right and left pointers of a node       *
*        and then free the current node using the free() C utility.        *
*        Note that all nodes attached to the node via right or left        *
*        pointers area also freed, so be sure that this is the desired     *
*        result when calling this routine.                                 *

*        This routine is supplied to allow successive calls to the         *
*        decmpinit routine.  It will free up the memory allocated          *
*        by previous calls to the decmpinit routine.  The call to free     *
*        a previous huffman tree is:  total = free_node(tree,(uint32_t) 0);    *
*        This call must be done by the programmer application routine      *
*        and is not done by any of these routines.                         *
*_HIST   16-AUG-89  Kris Becker U.S.G.S  Flagstaff Original Version        */
{
	if (pnode == (NODE *) NULL) return(total_free);
	
	if (pnode->right != (NODE *) NULL)
		total_free = free_node(pnode->right,total_free);
	if (pnode->left != (NODE *) NULL)
		total_free = free_node(pnode->left,total_free);

	free((char *) pnode);
	return(total_free + 1);
}
