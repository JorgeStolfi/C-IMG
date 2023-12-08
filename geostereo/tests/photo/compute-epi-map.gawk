#! /bin/gawk -f
# Last edited on 2003-02-16 16:37:12 by stolfi

BEGIN {
  usage = ( \
    "compute-epi-map -f r3x3.gawk \\\n" \
    "  [ -v width=WD ] [ -v check=BOOL ] \\\n" \
    "  < image.grmap > image.epmap" \
  );
  abort = -1;
  if (width == "") { width = 512; }
  height = width*384/512; 
  if (check == "") { check = 1; }
  split("",Md); split("",Mi);
  nlin = 0;
}

(abort >= 0) { exit abort; }

/^ *([#]|$)/ { next; }

(NF == 3) {
  # The map {M} takes full-image {H,V} (pix) to grid coords {X,Y} (mm).
  if (nlin < 3) 
    { i = nlin; 
      Md[i,0] = $1; Md[i,1] = $2; Md[i,2] = $3;
    }
  else
    { i = nlin - 3; 
      Mi[i,0] = $1; Mi[i,1] = $2; Mi[i,2] = $3;
    }
  nlin++;
  next;
}

//{ data_error(("bad line \"" $0 "\"")); }

END {
  if (abort >= 0) { exit abort; }
  if (nlin != 6) { data_error(("wrong line count = " nlin)); }
  # printf "Md:\n"; r3x3_print(Md);
  # Coordinate systems ({W = width}, {H = height}):
  #   CtrXY {X,Y}:   {[-W/2 _ +W/2] × [-H/2 _ +H/2]} (pix), {Y} up.
  #   FullHV {H,V}:  {[0 _ 2048] × [0 _ 1536]} (pix), {V} down.
  #   GridXY {X,Y}:  {[0 _ 250] × [0 _ 170]} (mm), {Y} up.
  # The matrices {Md,Mi} take FullHV to GridXY.
  # The pgmtran tool needs a matrix that maps CtrXY to CtrXY.
  split("",Hd); split("",Hi); # Maps CtrXY to FullHV.
  build_pre_scale_map(Hd,Hi);
  split("",Qd); split("",Qi); # Composite map {H*M}.
  r3x3_mul(Hd,Md,Qd);
  r3x3_mul(Mi,Hi,Qi);
  split("",Nd); split("",Ni); # Maps GridXY to CtrXY.
  build_pos_scale_map(Nd,Ni);
  split("",Pd); split("",Pi); # Composite map {H*M*N}.
  r3x3_mul(Qd,Nd,Pd);
  r3x3_mul(Ni,Qi,Pi);
  # r3x3_norm(Pd, 1.0);
  # r3x3_norm(Pi, 1.0);
  r3x3_print(Pd);
  printf "\n";
  r3x3_print(Pi);
  if (check) { test_matrices(Pd,Pi); }
}

# These procedure assume that caller has initialized the matrix and
# vector arguments {split("",*)}.

function build_pre_scale_map(Sd,Si,    i,j,scale,dh,dv)
{
  # Sets {S} to the 3x3 homogeneous coordinate matrix that takes
  # CtrXY coords (actual-size center-relative {X,Y}, in pixels)
  # to FullHV coords (full-size top-left-relative {H,V}, in pixels).
  scale = 2048/width;
  dh = int((width - 250.0*scale)/2.0 + 0.5);
  dv = int((height - 170.0*scale)/2.0 + 0.5);
  for (i = 0; i < 3; i++) { for (j = 0; j < 3; j++) { Sd[i,j] = 0.0; } }
  Sd[0,0] = 1.0;
  Sd[0,1] = 2048/2.0;
  Sd[0,2] = 1536/2.0;
  Sd[1,1] = scale;
  Sd[2,2] = -scale;
  # printf "Sd:\n"; r3x3_print(Sd);
  r3x3_inv(Sd,Si);
  # printf "Si:\n"; r3x3_print(Si);
  # r3x3_norm(Sd, 1);
  # r3x3_norm(Si, 1);
}

function build_pos_scale_map(Sd,Si,    i,j,scale,dh,dv)
{
  # Sets {M} to the 3x3 homogeneous coordinate matrix that takes
  # the GridXY coords (grid coordinates {[0_25] × [0_17]}, in mm)
  # to CtrXY coords (actual-size center-relative {X,Y}, in pixels), filling most 
  # of the actual image.
  
  scale = width/270.0;
  dh = int((width - 250.0*scale)/2.0 + 0.5);
  dv = int((height - 170.0*scale)/2.0 + 0.5);
  for (i = 0; i < 3; i++) { for (j = 0; j < 3; j++) { Sd[i,j] = 0.0; } }
  Sd[0,0] = 1.0;
  Sd[0,1] = - scale*250.0/2.0;
  Sd[0,2] = - scale*170.0/2.0;
  Sd[1,1] = scale;
  Sd[2,2] = scale;
  # printf "Sd:\n"; r3x3_print(Sd);
  r3x3_inv(Sd,Si);
  # printf "Si:\n"; r3x3_print(Si);
  # r3x3_norm(Sd, 1);
  # r3x3_norm(Si, 1);
}

function test_matrices(Pd,Pi,   ih,iv,p,q,rh,rv)
{ 
  split("", p);
  split("", q);
  printf "\n" > "/dev/stderr";
  printf "testing direct map\n" > "/dev/stderr";
  for(ih = 0; ih < 2; ih++)
    { for(iv = 0; iv < 2; iv++)
        { p[0] = 1; p[1] = ih*width; p[2] = -iv*height; 
          r3x3_mul_row(p,Pd,q);
          rh = q[1]/q[0];
          rv = q[2]/q[0];
          printf "  ( %7.1f %7.1f )  -> [ %7.2f %7.2f %7.2f ] = ( %7.2f %7.2f )\n", 
            p[1], p[2], q[0], q[1], q[2], rh, rv > "/dev/stderr";
        }
    }
  printf "\n" > "/dev/stderr";
  printf "testing inverse map\n" > "/dev/stderr";
  for(ih = 0; ih < 2; ih++)
    { for(iv = 0; iv < 2; iv++)
        { q[0] = 1; q[1] = ih*width; q[2] = -iv*height;
          r3x3_mul_row(q,Pi,p);
          rh = p[1]/p[0];
          rv = p[2]/p[0];
          printf " ( %7.1f %7.1f )  -> [ %7.2f %7.2f %7.2f ] = ( %7.2f %7.2f )\n", 
            q[1], q[2], p[0], p[1], p[2], rh, rv > "/dev/stderr";
        }
    }
  printf "\n" > "/dev/stderr";
}

function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}

function arg_error(msg)
{ printf "**%s\n", msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
