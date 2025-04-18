indent: Standard input:181: Warning:old style assignment ambiguity in "=*". Assuming "= *"

/*         Table of content        
FFT
FFT2d
Lanco
LiveAgua 
Pgmer
ReadM
Show 
SaveCT
SaveM
SaveT */

#include <math.h>
#include <stdio.h>
#include <malloc.h>
#include <strings.h>
typedef struct
{
  int m;
  int n;
  float *a;
} pgmDscr;
typedef struct
{
  int m;
  int n;
  int *num;
} MDscr;

/* reprint is used by SaveCT */

void
reprint (char *name, int length)
{
  char *rname;
  int i, jm, jn, m, n, l, k, trig;
  char cz;
  int ileft[64], iright[64];
  char *left[64], *right[64];
  FILE *f1, *f2;

/****** Openning files, allocating memory. ******/
  rname = calloc (32, sizeof (char));
  strcat (rname, name);
  strcat (rname, ".CAdat");
  f1 = fopen (rname, "r");
  free (rname);
  rname = calloc (32, sizeof (char));
  strcat (rname, name);
  strcat (rname, "CT.dat");
  f2 = fopen (rname, "w");
  fprintf (f2, "\n\n   *** THIS IS %s ***\n", name);
  free (rname);

/****** Reading dimention of the array ******/
  fscanf (f1, "%i %i\n", &m, &n);

  for (i = 0; i < n * 2; i++)
    {
      left[i] = calloc (8, sizeof (char));
      right[i] = calloc (8, sizeof (char));
    }

/****** Main loop ******/
  for (jm = 0; jm < m; jm++)	/* Loop N0 */
    {
      fprintf (f2, "\n");
      for (i = 0; i < n * (length + 2) + 3; i++)
	fprintf (f2, "-");
      fprintf (f2, "\n%.2i>", jm);

/* Loop N 1 */
      for (jn = 0; jn < n * 2; jn++)
	{
	  trig = 0;
	  fscanf (f1, "%c", &cz);

/****** Inner loop, N1, until CR is reached ******/
	  ileft[jn] = 0;
	  iright[jn] = 0;
	  for (i = 0; cz != '\n'; i++)
	    {

	      if (trig == 0)
		{
		  *(left[jn] + ileft[jn]) = cz;
		  ileft[jn]++;
		}
	      else
		{
		  *(right[jn] + iright[jn]) = cz;
		  iright[jn]++;
		};
	      if (cz == '.')
		trig = 1;
	      fscanf (f1, "%c", &cz);
	    }
/* End of inner loop */

	  *(right[jn] + iright[jn]) = '\0';
	  *(left[jn] + ileft[jn]) = '\0';

	}
/* End of N1 loop */
/* Loop N2 */
      for (jn = 0; jn < 2 * n; jn = jn + 2)
	{

	  fprintf (f2, "| ");

/****** Printing digits before . & dot itself ******/
	  for (k = 0; k < ileft[jn]; k++)
	    fprintf (f2, "%c", *(left[jn] + k));


/****** length-k digits left to print ******/

/****** Printing digits left ******/
	  for (i = 0; i < length - k & *(right[jn] + i) != '\0'; i++)
	    fprintf (f2, "%c", *(right[jn] + i));

/****** Putting spaces if neccessary ******/
	  for (l = 0; l < length - k - i; l++)
	    fprintf (f2, " ");

	}
/* End of loop N 2 */

      fprintf (f2, "| ");

      fprintf (f2, "\n   ");
/* Loop N3 */
      for (jn = 1; jn < 2 * n; jn = jn + 2)
	{

	  fprintf (f2, "| ");

/****** Printing digits before . & dot itself ******/
	  for (k = 0; k < ileft[jn]; k++)
	    fprintf (f2, "%c", *(left[jn] + k));

/****** length-k-1 digits left to print ******/

/****** Printing digits left ******/
	  for (i = 0; i < length - k & *(right[jn] + i) != '\0'; i++)
	    fprintf (f2, "%c", *(right[jn] + i));

/****** Putting spaces if neccessary ******/
	  for (l = 0; l < length - k - i; l++)
	    fprintf (f2, " ");

	}
/* End of Loop N 3 */

      fprintf (f2, "| ");

    }
/****** End of main loop ******/

  fprintf (f2, "\n");
  for (i = 0; i < n * (length + 2) + 3; i++)
    fprintf (f2, "-");

  fclose (f1);
  fclose (f2);
  free (left);
  free (right);

}

/* Records plr array to 2 pgm-files : ..A & ..F - ampl & phase,   */
/* dividing every member by factor maxampl. Phase range +-pi/2    */

void
SShow (char *name, int m, int n, plr * ax, float maxampl)
{
  FILE *f1, *f2;
  char *nname;
  plr z;
  int i, j;
  float trace1;
  int trace2;


  nname = calloc (20, sizeof (char));

  strcat (nname, name);
  strcat (nname, "A.pgm");
  f1 = fopen (nname, "w");
  free (nname);

  nname = calloc (20, sizeof (char));
  strcat (nname, name);
  strcat (nname, "F.pgm");
  f2 = fopen (nname, "w");

  fprintf (f1, "P2 %i %i 128", m, n);
  fprintf (f2, "P2 %i %i 128", m, n);
  free (nname);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      {
	z = *(ax + i * n + j);
	trace1 = z.r * 126.0 / maxampl;
	trace2 = (int) trace1;
	fprintf (f1, "\n%i", trace2);
	trace1 = 63.0 + z.an * 40.0;
	trace2 = (int) indent: Standard input:385: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:386: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:387: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:388: Warning:old style assignment ambiguity in "=*". Assuming "= *"

trace1;
	fprintf (f2, "\n%i", trace2);
      }

  fclose (f1);
  fclose (f2);
}


/****** Records cmp array to 2 pgm-files ******************************/

void
Show (char *name, int m, int n, cmp * x)
{
  plr *y;
  float mx;
/* cmp *extd;
int mm,nn,llm,lln;

llm=log2i(m); lln=log2i(n); */

  y = (plr *) calloc (m * n, sizeof (plr));

  mx = mCmpdd (m, n, x, y);

  SShow (name, m, n, y, mx);
}



/***** Saves matrix to dat-file ***********************************************/

void
SaveT (char *name, int m, int n, int *x)
{
  int i, j;
  FILE *f;
  char *rname;

  rname = calloc (20, sizeof (char));
  strcat (rname, name);
  strcat (rname, "T.dat");
  f = fopen (rname, "w");

  for (i = 0; i < m; i++)
    {
      fprintf (f, "\n");
      for (j = 0; j < n; j++)
	fprintf (f, " %.2i", *(x + i * n + j));
    };

  fclose (f);
  free (rname);

}

/* SaveM saves int matrix to M.dat file */

void
SaveM (char *name, int m, int n, int *x)
{
  int i, j;
  FILE *f;
  char *rname;

  rname = calloc (20, sizeof (char));
  strcat (rname, name);
  strcat (rname, "M.dat");
  f = fopen (rname, "w");

  fprintf (f, "%i %i", m, n);

  for (i = 0; i < m * n; i++)
    fprintf (f, "\n%.3i", *(x + i * n + j));

  fclose (f);
  free (rname);

}

/***** Reads pgm-file and returns pgm-descriptor ******************************/

pgmDscr
Readpgm (char *name)
{
  int i, j;
  int m, n, mx, z;
  char c;
  FILE *f;
  char *rname;
  pgmDscr res;

  rname = calloc (20, sizeof (char));
  strcat (rname, name);
  strcat (rname, ".pgm");
  f = fopen (rname, "r");
  free (rname);

  for (i = 0; i < 3; i++)
    fscanf (f, "%c", &c);

  fscanf (f, "%i %i %i", &n, &m, &mx);
  res.a = (float *) calloc (m * n, sizeof (float));
  res.m = m;
  res.n = n;
  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      {
	fscanf (f, "%i", &z);
	*(res.a + i * n + j) = 1.0 * z / mx;
      }

  printf ("\n pgmRead message... Tudo bem.");
  return (res);
}


/***** Makes a pbm-file from array x ****************************************/

void
pbmer (char *name, int m, int n, int *x)
{
  int i, j;
  char *nname;
  FILE *f;

  nname = calloc (20, sizeof (char));
  strcat (nname, name);
  strcat (nname, ".pbm");
  f = fopen (nname, "w");

  fprintf (f, "P1 %i %i", n, m);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      fprintf (f, "\n%i", *(x + i * n + j));

  fclose (f);
  free (nname);

}



/* ReadM reads matrix from M.dat-file and returns M-descriptor */

MDscr
ReadM (char *name)
{
  char *rname;
  FILE *f;
  MDscr res;
  int i, m, n;

  rname = calloc (20, sizeof (char));
  strcat (rname, name);
  strcat (rname, "M.dat");
  f = fopen (rname, "r");
  free (rname);

  fscanf (f, "%i %i", &m, &n);
  res.m = m;
  res.n = n;

  res.num = (int *) calloc (m * n, sizeof (int));

  for (i = 0; i < m * n; i++)
    fscanf (f, "%i", (res.num + i));

  fclose (f);
  printf ("\n ReadM message... Tudo bem. %i %i %i *", res.m, res.n,
	  *(res.num));
  return (res);
}

/***** Do not use it *******************************************************/

void
Filter (int m, int n, cmp * ax, cmp * ay, int par)
{
  int i, j, krit1, krit2;
  cmpindent: Standard input:409: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:432: Warning:old style assignment ambiguity in "=-". Assuming "= -"

indent: Standard input:437: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:480: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:567: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:568: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:578: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:588: Warning:old style assignment ambiguity in "=*". Assuming "= *"

indent: Standard input:601: Warning:old style assignment ambiguity in "=*". Assuming "= *"

 zero;

  zero.re = 0.0;
  zero.im = 0.0;
  printf ("\n Par= %i Zero: re %f im %f", par, zero.re, zero.im);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      {
	krit1 = n * (i - m / 2) + m * (j - n / 2);
	krit2 = n * (i - m / 2) - m * (j - n / 2);
	if ((krit1 > par) | (krit1 < -par) | (krit2 > par) | (krit2 < -par))
	  *(ay + i * n + j) = zero;
	else
	  *(ay + i * n + j) = *(ax + i * n + j);
      }

}


void
ReadMaskInt (char *name, int *ay)
{
  int i, j, m, n;
  int z;
  FILE *f;

  f = fopen (name, "r");

  fscanf (f, "%i %i", &m, &n);

  for (i = 0; i < m; i++)
    for (j = 0; j < n; j++)
      {
	fscanf (f, "%i", &z);
	*(ay + i * n + j) = z;
      }

  fclose (f);
}
