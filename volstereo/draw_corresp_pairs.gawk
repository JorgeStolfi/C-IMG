#! /usr/bin/gawk -f
# Last edited on 2003-04-17 15:58:14 by stolfi

BEGIN {
  abort = -1;
  if (frow == "") { arg_error("must specify frow"); }
  if (fcol == "") { arg_error("must specify fcol"); }
  if (rows == "") { arg_error("must specify rows"); }
  if (cols == "") { arg_error("must specify cols"); }
  if (scale == "") { scale = 1.0; }
  if (bothsides == "") { bothsides = 0; }
  npairs = 0;
}

(abort >= 0) { exit abort; }

/^ *([#]|$)/ { next; }

/[0-9]/ {
  if (NF != 3) { data_error("wrong corresp point format"); }
  col = $1; row1 = $2; row2 = $3;
  if ((col < fcol) || (col >= fcol + cols)) { next; }
  if ((row1 < frow) || (row1 >= frow + rows)) { next; }
  if ((row2 < frow) || (row2 >= frow + rows)) { next; }
  x1 = int(scale*(row1-frow)); 
  x2 = int(scale*(row2-frow));
  printf "circle %d,%d %d,%d\n", x1, x2, x1+3, x2;
  if (bothsides) { printf "circle %d,%d %d,%d\n", x2, x1, x2+3, x1; }
  npairs++;
  next;
}

END {
  if (abort >= 0) { exit abort; }
  # Ensure at least one command in output file ("convert" breaks otherwise):
  if (npairs == 0) 
    { mrow = int(frow + rows/2);
      printf "circle %d,%d %d,%d\n", mrow, mrow, mrow+1, mrow;
    }
}

function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}

function arg_error(msg)
{ printf "**%s\n", msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
