#! /bin/gawk -f
# Last edited on 2003-01-27 00:48:49 by stolfi

BEGIN {
  usage = ( \
    "compute-grid-map -f r3x3.gawk \\\n" \
    "  [ -v check=BOOL ] \\\n" \
    "  < image.pts > image.grmap" \
  );
  abort = -1;
  np = 0;
  if (check == "") { check = 1; }
  split("", x);
  split("", y);
  split("", z);
  split("", h);
  split("", v);
}

(abort >= 0) { exit abort; }

($6 ~ /^grid[0-9][0-9]/) {
  if (($5 + 0) != 0) { data_error(("bad Z = \"" $5 "\"")); }
  h[np] = $1;
  v[np] = $2;
  x[np] = $3;
  y[np] = $4;
  z[np] = $5;
  np++;
}

END {
  if (abort >= 0) { exit abort; }
  split("",Md); split("",Mi);
  split("",Nd); split("",Ni);
  split("",Pd); split("",Pi);
  find_proj_map(h,v,np,Md,Mi);
  find_proj_map(x,y,np,Nd,Ni);
  r3x3_mul(Mi,Nd,Pd);
  r3x3_mul(Ni,Md,Pi);
  r3x3_norm(Pd, 1.0);
  r3x3_norm(Pi, 1.0);
  r3x3_print(Pd);
  printf "\n";
  r3x3_print(Pi);
  
  # Testing the map:
  if (check) { test_proj_map(Pd,Pi,h,v,x,y,np); }
}

function test_proj_map(Pd,Pi,h,v,x,y,np,   i,j,k,p,q,rx,ry,rh,rv)
{ split("", p);
  printf "\n" > "/dev/stderr";
  printf "testing map (H,V)->(X,Y)\n" > "/dev/stderr";
  for(k = 0; k < np; k++)
    { p[0] = 1; p[1] = h[k]; p[2] = v[k]; 
      r3x3_mul_row(p,Pd,q);
      rx = q[1]/q[0];
      ry = q[2]/q[0];
      printf "  %6.1f %6.1f  -> %7.2f %7.2f   error  %7.3f %7.3f\n", 
        h[k], v[k], rx, ry, rx-x[k], ry-y[k] > "/dev/stderr";
    }
  printf "\n" > "/dev/stderr";
  printf "testing map (X,Y)->(H,V)\n" > "/dev/stderr";
  for(k = 0; k < np; k++)
    { q[0] = 1; q[1] = x[k]; q[2] = y[k]; 
      r3x3_mul_row(q,Pi,p);
      rh = p[1]/p[0];
      rv = p[2]/p[0];
      printf "  %6.1f %6.1f  -> %7.2f %7.2f   error  %7.3f %7.3f\n", 
        x[k], y[k], rh, rv, rh-h[k], rv-v[k] > "/dev/stderr";
    }
  printf "\n" > "/dev/stderr";
}

# These procedure assume that caller has initialized the matrix and
# vector arguments {split("",*)}.

function find_proj_map(X,Y,np,Md,Mi,    i,j,u,v,W)
{
  # Sets {M} to the 3x3 homogeneous coordinate matrix that takes
  # the canonical frame to the first four points {(X[i],Y[i])}. 
  
  if (np < 4) { data_error(("only " np " grid points --- needs at least 4")); }
  W = 1.0;
  for (i = 0; i < 3; i++) 
    { Md[i,0] = W; Md[i,1] = W*X[i]; Md[i,2] = W*Y[i]; }
  r3x3_inv(Md,Mi);
  # printf "Md:\n"; r3x3_print(Md);
  # printf "Mi:\n"; r3x3_print(Mi);
  split("",u);
  split("",v);
  u[0] = W; u[1] = W*X[3]; u[2] = W*Y[3];
  r3x3_mul_row(u,Mi,v);
  for (i = 0; i < 3; i++) 
    { for (j = 0; j < 3; j++)
        { Md[i,j] *= v[i]; Mi[j,i] /= v[i]; }
    }
  r3x3_norm(Md, 1);
  r3x3_norm(Mi, 1);
}

function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}

function arg_error(msg)
{ printf "**%s\n", msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
