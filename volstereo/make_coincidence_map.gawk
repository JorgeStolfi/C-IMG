#! /usr/bin/gawk -f
# Last edited on 2003-04-09 21:57:43 by stolfi

BEGIN {
  abort = -1;
  split("", smp);
  split("", dst);
  maxval_ot = 256*128-1;
  eps = 0.02;  # Assumed noise in [0_1] normalized intensities
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
      row = -1; col = width;
    }
  next;
}

/^ *[0-9]+/ {
  # Store samples:
  for (i = 1; i <= NF; i++)
    { s = $(i);
      if (s !~ /^[0-9]+/) { data_error("bad sample format"); }
      if (col >= samples_per_row) { row++; col = 0; }
      smp[row,col] = (s + 0.0)/maxval_in;
      col++;
    }
  next;
}
  
/./ { data_error(("bad line format \"" $0 "\"")); }

END {
  if (abort >= 0) { exit abort; }
  n = row;
  printf "P2\n"; 
  printf "%d %d\n", n, n; 
  printf "%d\n", maxval_ot; 
  eps2 = eps*eps;
  max_d = log((1 + eps2)/eps2);
  max_sum = samples_per_row*(max_d*max_d);
  for (i = 0; i < n; i++)
    { linsz = 0;
      for (j = 0; j < n; j++)
        { sum = 0.0; 
          for (col = 0; col < samples_per_row; col++)
            { a = smp[i,col]; b = smp[j,col];
              d = log((a*a + eps2)/(b*b + eps2));
              sum += d*d;
            }
          y = sqrt(sum/max_sum);
          yq = int(y*maxval_ot + 0.5);
          if (yq > maxval_ot) { yq = maxval_ot; }
          if (linsz >= 20) { printf "\n"; linsz = 0; }
          printf " %5d", yq;
          linsz++;
        }
      printf "\n";
    }
}

function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
