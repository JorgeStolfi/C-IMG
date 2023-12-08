// Last edited on 2003-04-17 00:48:43 by stolfi

#declare monolith_size = < +020.0, +010.0, +037.0 >;

// Monolith corners in its own coordinate system
#declare monolith_min = < -monolith_size.x/2, -monolith_size.y/2, 0.0 >; 
#declare monolith_max = monolith_min + monolith_size; 

#declare tx_monolith_side =
  texture {
    pigment { color rgb < 0.000, 0.550, 0.200 > }
    finish { ambient 0.0 diffuse 0.98 specular 0.1 roughness 0.1 }
  }

#declare tx_monolith_top =
  texture {
    pigment { color rgb < 0.980, 0.980, 0.980 > }
    finish { ambient 0.0 diffuse 1.0 }
  }

#declare monolith = 
  intersection {
    box {
      monolith_min,
      monolith_max + eps
      texture { tx_monolith_side }
    }
    plane { z, 0 
      texture { tx_monolith_top }
      translate monolith_max.z * z
    }
  }
