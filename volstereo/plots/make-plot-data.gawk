#! /usr/bin/gawk -f
# Last edited on 2003-04-10 00:42:46 by stolfi

BEGIN {
  # Reformats a single-column PPM/PGM (ascii) image into a data file
  # suitable for plotting by "gnuplot"
  abort = -1;
  hdstage = 0; # Header-parsing state
}

(abort >= 0) { exit abort; }

/^[#]/ { next; }

(hdstage < 4) { 
  i = 1;
  while ((hdstage < 4) && (i <= NF)) 
    { f = $(i); i++;
      if (hdstage == 0)
        { if (f !~ /^P[23]$/) { data_error("bad header"); }
          samples_per_pixel = (f == "P2" ? 1 : 3);
        }
      else if (hdstage == 1)
        { if (f !~ /^[0-9]+$/) { data_error("bad header"); }
          width = f + 0;
          
        }
      else if (hdstage == 2)
        { if (f !~ /^[0-9]+$/) { data_error("bad header"); }
          height = f + 0;
        }
      else if (hdstage == 3)
        { if (f !~ /^[0-9]+$/) { data_error("bad header"); }
          maxval_in = f + 0;
        }
      else
        { data_error("bad hdstage"); }
      hdstage++;
    }
  if (hdstage == 4)
    { if (i <= NF) { data_error("garbage in header"); }
      samples_per_row = samples_per_pixel * width;
      row = 0; col = 0;
    }
  next;
}

/^ *[0-9]+/ {
  # Print samples
  for (i = 1; i <= NF; i++)
    { s = $(i);
      if (s !~ /^[0-9]+/) { data_error("bad sample format"); }
      if (col >= samples_per_row) { printf "\n"; row++; col = 0; }
      printf " %5d", s;
      col++;
    }
  next;
}
  
/./ { data_error(("bad line format \"" $0 "\"")); }

END {
  if (abort >= 0) { exit abort; }
  if (col != samples_per_row) { data_error("unexpected EOF"); }
  printf "\n"; 
}

function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
