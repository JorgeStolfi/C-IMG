#! /usr/bin/gawk -f
# Last edited on 2023-02-08 09:31:28 by stolfi
# To be included in parms-to-geomview and ppmhist-yo-geomview 

function gamma_in(y, gm)
{
  if (y == 0) 
    { return y; }
  else if (y == 1)
    { return y; }
  else if (gm == 1) 
    { return y; }
  else
    { return exp(log(y)/gm); }
}

function log_scale_in(y,eps,   lmag,yabs,ylog)
{
  if (y != 0.0)
    { lmag = (eps < 1.0 ? 1.0/eps/(1 - log(eps)) : 1.0);
      yabs = (y < 0.0 ? -y : y);
      if (yabs > eps)
        { ylog = eps*(log(yabs/eps) + 1); 
          y = (y < 0.0 ? -ylog : ylog);
        }
      y *= lmag;
    }
  return y;
}

