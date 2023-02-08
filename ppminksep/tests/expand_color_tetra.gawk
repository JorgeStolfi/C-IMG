#! /usr/bin/gawk -f
# Last edited on 2004-07-29 19:48:48 by stolfi

BEGIN {
  usage = ( ARGV[0] " -v factor=NUM < IN.parms > OUT.parms" ); 
  
  # Reads a parmsfile, outputs the same with the color tetrahedron
  # expanded by the given factor.
  if (factor == "") { arg_error(("must define \"factor\"")); }
  
  split("", fname); 
  split("", r); split("", g); split("", b);
  split("", cmt);
  n = 0;
}

/^[ ]*[-]layer/ {
  if (match($0, /[ ]*[\#]/))
    { cmt[n] = substr($0,RSTART);
      $0 = substr($0, 1, RSTART-1); 
    }
  else
    { cmt[n] = 0; }
  fname[n] = $2;
  if ($3 != "-color") { data_error(("bad -layer fmt")); }
  if (NF == 6)
    { den = 1; }
  else if (NF == 8)
    { if ($7 != "/") { data_error(("bad -layer fmt")); }
      den = $8 + 0;
    }
  else
    { data_error(("bad -layer fmt")); }
  r[n] = $4/den;
  g[n] = $5/den;
  b[n] = $6/den;
  n++;
  if (n == 4)
    { sr = sg = sb = 0;
      for (i = 0; i < n; i++)
        { sr += r[i]; sg += g[i]; sb += b[i]; }
      sr /= n; sg /= n; sb /= n;
      den = 255;
      for (i = 0; i < n; i++)
        { r[i] = int(den*(sr + factor*(r[i] - sr)) + 0.5);
          g[i] = int(den*(sg + factor*(g[i] - sg)) + 0.5);
          b[i] = int(den*(sb + factor*(b[i] - sb)) + 0.5);
          if (r[i] < 0) { r[i] = 0; } if (r[i] > den) { r[i] = den; }
          if (g[i] < 0) { g[i] = 0; } if (g[i] > den) { g[i] = den; }
          if (b[i] < 0) { b[i] = 0; } if (b[i] > den) { b[i] = den; }
          printf "-layer %s -color %03d %03d %03d / %03d%s\n", \
            fname[i], r[i], g[i], b[i], den, cmt[i];
        }
    }
  next;
}

// { print; }

END {
  if (n != 4) { data_error(("wrong num of layers")); }
}

function data_error(msg)
{ 
  printf "line %d: %s\n", FNR, msg > "/dev/stderr";
  abort = 1;
  exit 1
}

      
      
          
