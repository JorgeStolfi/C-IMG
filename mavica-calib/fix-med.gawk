#! /usr/bin/gawk -f
# Last edited on 2003-11-20 00:45:33 by stolfi

BEGIN {
  split("", fld);
}

/^ *[0-9]+[ ]/ {
  if ($1 < 0)
    { printf "bad rec" > "/dev/stderr"; }
  else if ($1 < 20)
    { i = $1;
      for (j = 0; j < 4; j++) { fld[i,j] = $(j+2); }
    }
  else if ($1 < 40)
    { i = $1 - 20;
      for (j = 0; j < 3; j++) { fld[i,4+j] = $(j+2); }
    }
  else
    { printf "bad rec" > "/dev/stderr"; }
}

END {
  for (i = 0; i < 20; i++)
    { printf "%02d %4.2f", i, fld[i,0];
      printf "  %6.2f %6.2f %6.2f", fld[i,1], fld[i,2], fld[i,3];
      printf "  %6.2f %6.2f %6.2f", fld[i,4], fld[i,5], fld[i,6];
      printf "\n";
   }
}
