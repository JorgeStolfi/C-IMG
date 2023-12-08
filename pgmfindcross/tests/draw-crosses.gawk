#! /bin/gawk -f
# Last edited on 2003-04-07 03:58:59 by stolfi

BEGIN {
  usage = ( \
    "draw-crosses \\\n" \
    "  [ -v scale=NUM ] \\\n" \
    "  [ -v radius=NUM ] \\\n" \
    "  < image.crs > image.drw" \
  );
  
  # Reads a list of line crossings, one per line, each given by its
  # {H} and {V} coordinates and the line directions {R} and {S}.
  # Outputs a string of plotting commands, suitable as the "-draw"
  # option of "convert", that draws crosshais at specified points.
  # Assumes that the commands will be applied to a copy of the image
  # that has been scaled by {scale}.

  if (scale == "") { scale = 1.0; }
  if (radius == "") { radius = 5.0; }
  abort = -1;
}

(abort >= 0) { exit abort; }

/^ *[-+]?[0-9]+([.][0-9]*)? / {
  if (NF != 4) { data_error("invalid format"); }
  h = $1;
  v = $2;
  a = $3;
  b = $4;
  plot_point(h,v);
  plot_line(h,v,a + 0.0);
  plot_line(h,v,a + 180.0);
  plot_line(h,v,b + 0.0);
  plot_line(h,v,b + 180.0);
}

function plot_point(h,v,  hc,vc,hr,vr)
{
  # Plots a dot at {h,v}.
  hc = int(scale*h);
  hr = hc + int(radius + 0.5);
  vc = int(scale*v);
  vr = vc;
  printf "circle %d,%d %d,%d\n", hc,vc, hr,vr;
}

function plot_line(h,v,a,  ar,ca,sa,h1,v1,h2,v2)
{
  # Plots a line segment out of {h,v} in the direction {a}.
  ar = 3.1415926*a/180.0;
  ca = cos(ar); sa = sin(ar);
  printf "a = %.3f  ca = %.8f sa = %.8f\n", a, ca, sa > "/dev/stderr";
  
  h1 = int(scale*h + 2*radius*ca); 
  v1 = int(scale*v + 2*radius*sa);
  
  h2 = int(scale*h + 4*radius*ca); 
  v2 = int(scale*v + 4*radius*sa);
  printf "line %d,%d %d,%d\n", h1,v1, h2,v2;
}

END {
  if (abort >= 0) { exit abort; }
}

function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}

function arg_error(msg)
{ printf "**%s\n", msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
