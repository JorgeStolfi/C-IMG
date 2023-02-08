#! /usr/bin/gawk -f
# Last edited on 2004-07-30 15:59:02 by stolfi

BEGIN {
  # Reads a ".parms" file, turns it into a geomview model showing the four colors
  usage = ( ARGV[0] " -f color-conv.gawk\\\n" \
    "  [ -v gmr=N ] [ -v gmg=N ] [ -v gmb=N ] \\\n" \
    "  [-v logScale=BOOL] \\\n" \
    "  < file.parms > file.gvtop" \
  );

  if (gmr == "") { gmr = 1.0; }
  if (gmg == "") { gmg = 1.0; }
  if (gmb == "") { gmb = 1.0; }
  if (logScale == "") { logScale = 0; }

  printf "{\n"
  printf "  LIST\n"
  k = 0;
  split("", rt); split("", gt); split("", bt); # Layer colors
  split("", rm); split("", gm); split("", bm); # Stretch colors
  eps = 1.0/255.0; # Log scale threshhold
  stretched = 0; # True if stretched image was requested.
  map_given = 1; # For now; true if "-map" options present.
}

# Get rid of comments
/[\#]/ { gsub(/[\#].*$/, "", $0); }

# Should we check whether the -inGamma and -logScale options in the
# ".parms" file are consistent with those given on the command line?
# 
# /^ *[-]inGamma /{
#   # Input gamma line
#   gmr = $2; gmg = $3; gmb = $4;
#   printf "using gamma = %6.4f %6.4f %6.4f\n", gmr, gmg, gmb > "/dev/stderr";
#   next;
# }
# 
# /^ *[-]logScale /{
#   # Log scale option
#   logScale = 1;
#   printf "using log scale\n" > "/dev/stderr";
#   next;
# }

/^ *[-]stretched/ {
  stretched = 1;
  next;
}
  
/^ *[-]layer / {
  # Layer line
  icolor = 3;
  if ($(icolor) != "-color")
    { printf "bad layer line [%s]\n", $0 > "/dev/stderr"; exit 1; }
  # Get "-color" information:
  imap = icolor+1 + grab_color(icolor+1, rt, gt, bt, k);

  # Get "-map" information: 
  if ((imap <= NF) && ($(imap) == "-map"))
    { iend = imap+1 + grab_color(imap+1, rm, gm, bm, k); }
  else
    { map_given = 0; }
  
  # Check format:
  if (NF != iend-1)
    { printf "bad layer NF [%s]\n", $0 > "/dev/stderr"; exit 1; }
 
  k++;
  next;
}

END {
  # Show tetrahedron edges
  if (k != 4) 
    { printf "found %d layers - aborted\n", k > "/dev/stderr"; exit 1; }
  show_tetrahedron(rt,gt,bt, 0.2,0.5,0.2);
  if (stretched && map_given)
    { show_tetrahedron(rm,gm,bm, 0.7,0.3,0.7); }
  show_rgb_cube();
  printf "}\n"; # End of LIST
  printf "< bola.off\n";
}

function grab_color(i,rv,gv,bv,k,  r,g,b,nf,den)
{
  # Parses fields {$(i)} through {$(i+2)} or {$(i+4)} as an RGB spec,
  # maps it though the current gamma and logscale options,
  # saves it in {rv[k],gv[k],bv[k]}.
  # Returns number of fields parsed
  
  if ($(i+3) == "/")
    { den = $(i+4) + 0.0; nf = 5; }
  else
    { den = 1.0; nf = 3; }
  if (i+nf-1 > NF) 
    { printf "bad layer NF [%s]\n", $0 > "/dev/stderr"; exit 1; }
  r = gamma_in($(i+0)/den, gmr); 
  g = gamma_in($(i+1)/den, gmg);
  b = gamma_in($(i+2)/den, gmb);
  if (logScale) 
    { r = log_scale_in(r,eps); 
      g = log_scale_in(g,eps); 
      b = log_scale_in(b,eps); 
    }
  rv[k] = r;
  gv[k] = g;
  bv[k] = b;
  
  return nf;
}

function show_color(r,g,b,size)
{
  x = r - 0.5;
  y = g - 0.5;
  z = b - 0.5;
  radius = size/255.0;
  printf "    { \n";
  printf "      appearance { material { ka 0.7 ambient %6.4f %6.4f %6.4f kd 0.3 diffuse %6.4f %6.4f %6.4f } }\n", r, g, b, r, g, b;
  printf "      INST geom { : bola }\n"; 
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

function show_tetrahedron(rv,gv,bv, re,ge,be,   k)
{
  for (k = 0; k < 4; k++)
    { show_color(rv[k],gv[k],bv[k], 5.0); }
  printf "    {\n"
  set_edge_appearance(re,ge,be);
  printf "      SKEL\n"
  printf "        4 6\n"
  for (k = 0; k < 4; k++)
    { printf "        %7.4f %7.4f %7.4f\n", rv[k]-0.5, gv[k]-0.5, bv[k]-0.5; }
  printf "        2 0 1\n";
  printf "        2 0 2\n";
  printf "        2 0 3\n";
  printf "        2 1 2\n";
  printf "        2 1 3\n";
  printf "        2 2 3\n";
  printf "    }\n";
}

function show_rgb_cube(   rc,gc,bc)
{
  for(rc = 0; rc <= 1; rc++)
    for(gc = 0; gc <= 1; gc++)
      for(bc = 0; bc <= 1; bc++)
        { show_color(rc,gc,bc, 3.0); }
  
  printf "    {\n"
  set_edge_appearance(0.25, 0.25, 0.25);
  printf "      SKEL\n"
  printf "        8 13\n"
  for(rc = 0; rc <= 1; rc++)
    for(gc = 0; gc <= 1; gc++)
      for(bc = 0; bc <= 1; bc++)
        { printf "        %7.4f %7.4f %7.4f\n", rc-0.5, gc-0.5, bc-0.5; }
  printf "        2 0 1\n";
  printf "        2 2 3\n";
  printf "        2 4 5\n";
  printf "        2 6 7\n";
  printf "        2 0 2\n";
  printf "        2 1 3\n";
  printf "        2 4 6\n";
  printf "        2 5 7\n";
  printf "        2 0 4\n";
  printf "        2 1 5\n";
  printf "        2 2 6\n";
  printf "        2 3 7\n";
  printf "        2 0 7\n";
  printf "    }\n";
}

function set_edge_appearance(r,g,b)
{
  printf "      appearance {\n";
  printf "        material {\n";
  printf "          ka 0.5 ambient 0.5 1.0 1.0\n";
  printf "          kd 0.5 diffuse 0.5 1.0 1.0\n";
  printf "          ks 0.0\n";
  printf "          edgecolor %6.4f %6.4f %6.4f \n", r, g, b;
  printf "        }\n";
  printf "      }\n"; 
}
