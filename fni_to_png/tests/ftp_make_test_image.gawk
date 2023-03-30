#! /usr/bin/gawk -f 
# Last edited on 2021-08-24 23:38:32 by stolfi

# Creates a FNI image suitable for testing {fni_to_png}. 

BEGIN { 
  nx = 64;
  ny = 48; 
  nc = 4; 
  rad = 20;
  
  printf "begin float_image_t (format of 2006-03-25)\n"; 
  printf "NC = %d\n", nc;
  printf "NX = %d\n", nx;
  printf "NY = %d\n", ny;

  for (iy = 0; iy < ny; iy++) { 
    for (ix = 0; ix < nx; ix++) { 
      printf "%5d %5d", ix, iy;
      for (ic = 0; ic < nc; ic++) { 
        ctrx = nx/2 + (2*ic/(nc-1)-1)*rad/4;
        ctry = ny/2;
        dx = ix + 0.5 - ctrx;
        dy = iy + 0.5 - ctry;
        d = sqrt(dx*dx + dy*dy);
        if (d > rad) 
          { v = -16; }
        else
          { v = (dx+dy)/rad*sqrt(rad*rad - d*d); }
        printf " %+11.7f", v;
      }
      printf "\n";
    }
    printf "\n";
  }
  printf "end float_image_t\n";
}
  
