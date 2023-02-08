#! /usr/bin/gawk -f
# Last edited on 2004-07-30 03:48:29 by stolfi

BEGIN {
  # Reads a "ppmhist" output, turns it into a geomview model of little spheres
  usage = ( ARGV[0] " -f color-conv.gawk\\\n" \
    "  [ -v gmr=N ] [ -v gmg=N ] [ -v gmb=N ] \\\n" \
    "  [-v logScale=BOOL] \\\n" \
    "  < file.parms > file.gvtop" \
  );
  
  if (gmr == "") { gmr = 1.0; }
  if (gmg == "") { gmg = 1.0; }
  if (gmb == "") { gmb = 1.0; }
  if (logScale == "") { logScale = 0; }
  printf "using gamma = %6.4f %6.4f %6.4f\n", gmr, gmg, gmb > "/dev/stderr";
  if (logScale) { printf "using log scale\n" > "/dev/stderr"; }
  eps = 1.0/255.0; # Log scale threshhold
  
  printf "{\n"
  printf "  LIST\n"
}

/^ *[0-9]/ {
  # Data line
  r = gamma_in($1/255.0, gmr);
  g = gamma_in($2/255.0, gmg);
  b = gamma_in($3/255.0, gmb);
  if (logScale) 
    { r = log_scale_in(r, eps); 
      g = log_scale_in(g, eps); 
      b = log_scale_in(b, eps); 
    }
  rawlum = $4/255.0; # Uncorrected luminosity from histogram 
  ct = $5;
  x = r - 0.5;
  y = g - 0.5;
  z = b - 0.5;
  # radius = (log(ct) + 1)/255.0;
  radius = (log(ct)/3 + 2)/255.0;
  printf "    { \n";
  printf "      appearance { material { ambient %6.4f %6.4f %6.4f } }\n", r, g, b;
  printf "      INST geom { : pingo }\n"; 
  printf "      "; print_transform(x,y,z,radius); printf "\n";
  printf "    }\n";
}

function print_transform(x,y,z,radius)
{
  printf "transform {";
  printf " %6.4f 0 0 0 ", radius;
  printf " 0 %6.4f 0 0 ", radius;
  printf " 0 0 %6.4f 0 ", radius;
  printf " %7.4f  %7.4f %7.4f 1 ", x, y, z;
  printf "}";
}

END {
  printf "}\n";
  printf "< pingo.off\n";
}
