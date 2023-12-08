#! /usr/bin/gawk -f
# Last edited on 2003-04-17 12:44:10 by stolfi

BEGIN {
  abort = -1;
  make_shape();
  rmax = 2.0;
  make_grains(rmax);
  write_model();
}

function make_shape()
{ 
  # Creates a blob that defines the shape of the fragment. Defines
  # {xblob,yblob,zblob,rblob,nblobs}, and also the size of the
  # bounding box {xsize,ysize,zsize} (starts at origin).
  
  split("", xblob);
  split("", yblob);
  split("", zblob);
  split("", rblob);
  # The blobs:
  xblob[0] = 05.00; yblob[0] = 12.00; zblob[0] = 08.00; rblob[0] = 05.00;
  xblob[1] = 15.00; yblob[1] = 08.00; zblob[1] = 06.00; rblob[1] = 10.00;
  xblob[2] = 12.00; yblob[2] = 09.00; zblob[2] = 04.00; rblob[2] = 09.00;
  nblobs = 3;
  # Maximum dimensions of fragment:
  xsize = 30;
  ysize = 27;
  zsize = 12;
}

function make_grains(rmax,   \
  rmin,step,jitt,n,xmin,xmax,xstep,ymin,ymax,ystep,zmin,zmax,zstep, \
  yshift,xshift,depth,darkfactor,zi,yi,xi,dx,dy,dz,x,y,z)
{
  # Generates a set of spherical grains in {xgr,ygr,zgr,rgr,Rgr,Ggr,Bgr}
  # indexed from 0 to {ngrains-1}.
  
  split("", xgr); # X of grain center
  split("", ygr); # Y of grain center
  split("", zgr); # Z of grain center
  split("", rgr); # Radius of grain
  split("", Rgr); # Hue of grain
  split("", Ggr); # Saturation of grain 
  split("", Bgr); # Value of grain 
  
  step = 1.20*rmax;       # Nominal step between grains
  
  # rmin = 0.75*rmax;       # Min grain radius
  # jitt = 2 * rmin - step; # Jitter radius for centers
  
  rmin = rmax;
  jitt = 0;
  
  n = 0;
  xmin = -xsize/2; xmax = xmin + xsize; xstep = step;
  ymin = -ysize/2; ymax = ymin + ysize; ystep = step * sqrt(3/4);
  zmin = rmax + jitt; zmax = zmin + zsize; zstep = step * sqrt(2/3);
  yshift = 0;
  depthnorm = zsize*zsize/4;
  for (zi = zmin; zi <= zmax; zi += zstep)
    { xshift = 0;
      depth = (zi - zmin)*(zmax - z)/depthnorm;
      darkfactor = 1.0 - 0.3*depth;
      for (yi = ymin; yi <= ymax; yi += ystep)
        { for (xi = xmin; xi <= xmax; xi += xstep)
            { do
                { dx = 2*rand()-1; dy = 2*rand()-1; dz = 2*rand()-1; }
              while (dx*dx + dy*dy + dz*dz > 1.0);
              x = xi + xshift + dx*jitt;
              y = yi + yshift + dy*jitt;
              z = zi + zshift + dz*jitt;
              if (inside(x,y,z))
                { xgr[n] = x;
                  ygr[n] = y;
                  zgr[n] = z;
                  rgr[n] = rmin + rand()*(rmax-rmin);
                  Rgr[n] = darkfactor*(0.80 + 0.15*rand());
                  Ggr[n] = darkfactor*(0.70 + 0.10*rand());
                  Bgr[n] = darkfactor*(0.30 + 0.30*rand());
                  n++;
                }
            }
          xshift = xstep/2 - xshift;
        }
      yshift = 2*ystep/3 - yshift;
    }
  ngrains = n;
}

function inside(x,y,z, fsum,i,dx,dy,dz,d2,r2)
{
  fsum = 0.0;
  for (i = 0; i < nblobs; i++) 
    { dx = x - xblob[i]; 
      dy = y - yblob[i];
      dz = z - zblob[i];
      d2 = dx*dx + dy*dy + dz*dz;
      r2 = rblob[i]*rblob[i];
      fsum += exp(-d2/r2);
    }
  return (fsum > 0.367879);
}

function write_model(  i)
{
  printf "  union {\n";
  for (i = 0; i < ngrains; i++)
    { printf "    sphere { < %6.2f, %6.2f, %6.2f >", xgr[i], ygr[i], zgr[i]; 
      printf ", %6.2f", rgr[i]; 
      printf " texture {"; 
      printf " pigment { color rgb <%4.2f, %4.2f, %4.2f > }", Rgr[i], Ggr[i], Bgr[i]; 
      printf " finish { ambient 0 diffuse 1 }"; 
      printf " }"; 
      printf " }\n"; 
    }
  printf "  }\n";
}

(abort >= 0) { exit abort; }


END {
  if (abort >= 0) { exit abort; }
}

function data_error(msg)
{ printf "%d: **%s\n", FNR, msg > "/dev/stderr"; 
  abort = 1; exit abort;
}
