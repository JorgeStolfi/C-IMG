// Last edited on 2003-06-21 16:51:15 by ra007998

#include "colors.inc"

background{ White }

#declare tx =
  texture {
    pigment { Black }
    finish { diffuse 0 ambient 1 }
  }

#declare camd = 3.6;
camera {
  location  1.0 * camd * < 0.0, 10.00, 0.0 >
  right     -0.50*x  
  up         0.50*y
  sky       z
  look_at   <  0.00, 0.00, 0.00 >
}
// cylinder { <0,0,0>, <0.5*camd,0,0>, 0.01*camd pigment{ Red }}
// cylinder { <0,0,0>, <0,0.5*camd,0>, 0.01*camd pigment{ Green }}
// cylinder { <0,0,0>, <0,0,0.5*camd>, 0.01*camd pigment{ Blue }}

#declare rmin = 0.03;

sphere { <0,0,0>, rmin 
  texture { pigment { color rgb <1,0,0> } finish { ambient 1 } }
}

cylinder { <-8,0,0>, <8,0,0>, rmin rotate 090*y texture { tx } }
