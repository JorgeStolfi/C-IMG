#! /usr/bin/gawk -f
# Last edited on 2003-04-15 13:04:53 by stolfi

BEGIN {
  abort = -1;
  if (img == "") { arg_error("must specify \"img\""); }
  if (col == "") { arg_error("must specify \"col\""); }
  if (irow0 == "") { arg_error("must specify \"irow0\""); }
  if (frow0 == "") { arg_error("must specify \"frow0\""); }
  if (irow1 == "") { arg_error("must specify \"irow1\""); }
  if (frow1 == "") { arg_error("must specify \"frow1\""); }
  if (channels == "") { arg_error("must specify \"channels\""); }
  
  rows0 = frow0 - irow0;
  rows1 = frow1 - irow1;
  rows = (rows0 < rows1 ? rows0 : rows1);
  printf "set terminal postscript eps color \"courier\" 12\n";
  printf "set output \"%s-%d.eps\"\n", img, col;
  printf "set size 2,1.5\n", img, col;
  printf "xmin = %d\n", -3 ;
  printf "xmax = %d\n", rows+2;
  printf "set xrange [(xmin):(xmax)]\n";
  sep = "plot";
  for (eye = 0; eye <= 1; eye++)
    { for (i = 1; i <= channels; i++)
        { printf "%s \\\n", sep; 
          sep = ",";
          printf "   \"%s-%s-%d.data\"", img, col, eye;
          xcoord = ( eye == 0 ? "column(0)" : "xmax - column(0)");
          yoffset = 100*(i-1);
          printf " using (%s):(column(%d)+%d)", xcoord, i, yoffset;
          printf " title \"eye %d", eye;
          if (channels > 1) { printf " chn %d", i; }
          printf "\"";
          printf " with linespoints lt %d", eye;
        }
    }
  printf "\n";
  printf "quit\n";
}

(abort >= 0) { exit abort; }


function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}

function arg_error(msg)
{ printf "**%s\n", msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
