// Last edited on 2010-04-29 21:04:43 by stolfilocal

global_settings { assumed_gamma 2.2 }

#include "colors.inc"
#include "shapes.inc"

background { color rgb <1.000, 1.000, 1.000>  }
#include "scene.inc"

#declare maxSize = max(maxX,maxY);

// Radius (half-diagonal) of domain rectangle:
#declare radiusXY = sqrt(maxX*maxX + maxY*maxY)/2;

// Desired approximate height of plot:
#declare plotZ = 0.50*radiusXY;

// Scale factor for Z coordinate:
#declare scaleZ = plotZ/scaleZ;

#declare color1_pos = < 1.000, 0.800, 0.700 >;
#declare color2_pos = < 1.000, 0.900, 0.800 >;

#declare gray1 = 0.3*color1_pos.x + 0.6*color1_pos.y + 0.1*color1_pos.z;
#declare gray2 = 0.3*color2_pos.x + 0.6*color2_pos.y + 0.1*color2_pos.z;

#declare color1_neg = gray1 - color1_pos;
#declare color2_neg = gray2 - color2_pos;

#declare tx_pos = 
  texture{ 
    pigment{ 
      checker color rgb color1_pos, color rgb color2_pos
      translate 0.5*y
      rotate 90*y
      scale < 100000, 1, stepZ*scaleZ >
    }  
    finish { diffuse 0.8 ambient 0.2 } 
  }

#declare tx_neg = 
  texture{ 
    pigment{ 
      checker color rgb color1_neg, color rgb color2_neg
      translate 0.5*y
      rotate 90*y
      scale < 100000, 1, stepZ >
    }  
    finish {diffuse 0.8 ambient 0.2} 
  }

#declare tx_zer = 
  texture{ 
    pigment{ color rgb < 0.750, 0.750, 0.750 > }
    finish {diffuse 0.8 ambient 0.2} 
  }

union{
  // Positive part:
  intersection{ 
    object{ scene_mesh scale < 1, 1, scaleZ > }
    plane{ -z, 0 translate +0.00001*plotZ }
    texture { tx_pos scale < 1, 1, scaleZ > } 
  }
  // Negative part:
  intersection{ 
    object{ scene_mesh scale < 1, 1, scaleZ > }
    plane{ +z, 0 translate -0.00001*plotZ }
    texture { tx_neg scale < 1, 1, scaleZ > } 
  }
  // Null part:
  intersection{ 
    object{ scene_mesh scale < 1, 1, scaleZ > }
    plane{ -z, 0 translate -0.00001001*plotZ }
    plane{ +z, 0 translate +0.00001001*plotZ }
    texture { tx_zer } 
  }

  translate <-(maxX)/2.0,-(maxY)/2.0,0 >
}

light_source { <20000, 50000, 100000> color rgb 0.80*<1,1,1> }
light_source { <-20000, 50000, 100000> color rgb 0.30*<1,1,1> }

#declare camera_dir = vrotate(vrotate(<1,0,0>, -viewElev*y), -viewAzim*z);
#declare camera_disp = viewDist*camera_dir;
#declare camera_radius = 0.95*radiusXY;
#declare camera_focus = <0,0,0>;
#declare image_rel_width = image_width;
#declare image_rel_height = image_height;

#include "camera.inc"
