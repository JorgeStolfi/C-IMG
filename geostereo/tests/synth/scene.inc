//  A test scene for stereo reconstruction
//  Last edited on 2003-03-05 10:58:47 by stolfi

// Dimensions in millimeters

background{ color rgb <0.00, 0.05, 0.30> }

light_source {
  <+750.000, +300.000, +600.000>
  color rgb <0.8, 0.8, 0.8>
} 

light_source {
  <-800.000, +500.000, +150.000>
  color rgb <0.4, 0.4, 0.4>
} 

light_source {
  <+250.000, +500.000, -400.000>
  color rgb <0.3, 0.3, 0.3>
} 

#declare color_rock = < 1.0, 0.95, 0.85 >;

#declare texture_rock = 
  texture{ 
    pigment{ color rgb color_rock }
    finish{ ambient 0.1 diffuse 0.9 }
    normal { granite 3.0 }
  }

union{
   torus{ 35.0, 12.0 }
   plane{ <0,1,0>, 0 }
   texture { texture_rock scale 5.0 }
 }

#declare EyeSep = 2.0;

#declare camDist = 30.0;
