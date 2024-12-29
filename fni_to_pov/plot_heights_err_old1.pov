// Last edited on 2010-04-29 22:38:21 by stolfilocal

global_settings { assumed_gamma 2.2 }

#include "colors.inc"
#include "shapes.inc"

background { color rgb <1.000, 1.000, 1.000>  }
#include "scene.inc"

// True to paint pos & neg with different shades.
#declare bicolor = false;

// Radius (half-diagonal) of domain rectangle:
#declare radiusXY = sqrt(maxX*maxX + maxY*maxY)/2;

#declare color1_pos = < 1.000, 0.500, 0.300 >;
#declare color2_pos = < 1.000, 0.900, 0.600 >;

#declare gray1 = 0.3*color1_pos.x + 0.6*color1_pos.y + 0.1*color1_pos.z;
#declare gray2 = 0.3*color2_pos.x + 0.6*color2_pos.y + 0.1*color2_pos.z;

#declare color1_neg = gray1 - color1_pos;
#declare color2_neg = gray2 - color2_pos;

#declare tx_scale = < 10*radiusXY, 10*radiusXY, stepZ >;

#declare clip_pos = 
  box{
    < -0.55*maxX, -0.55*maxY, 0.00001*radiusXY >,
    < +0.55*maxX, +0.55*maxY, 0.1*radiusXY + maxZ >
  }

#declare tx_pos = 
  texture{ 
    pigment{ 
      checker color rgb color1_pos, color rgb color2_pos
      translate < 0.5, 0.5, 0.0 >
      scale tx_scale
    }  
    finish { diffuse 0.8 ambient 0.2 } 
  }

#declare clip_neg = 
  box{
    < -0.55*maxX, -0.55*maxY, -0.1*radiusXY - maxZ >,
    < +0.55*maxX, +0.55*maxY, -0.00001*radiusXY >
  }

#declare tx_neg = 
  texture{ 
    pigment{ 
      checker color rgb color1_neg, color rgb color2_neg
      translate < 0.5, 0.5, 0.0 >
      scale tx_scale
    }  
    finish {diffuse 0.8 ambient 0.2} 
  }

#declare clip_zer = 
  box{
    < -0.55*maxX, -0.55*maxY, -0.00001*radiusXY >,
    < +0.55*maxX, +0.55*maxY, +0.00001*radiusXY >
  }

#declare tx_zer = 
  texture{ 
    pigment{ color rgb < 0.750, 0.750, 0.750 > }
    finish { diffuse 0.8 ambient 0.2 } 
  }

union{
  #if (bicolor)
    // Positive part:
    object{ 
      scene_mesh scale < 1, 1, scaleZ >
      clipped_by{ clip_pos }
      texture { tx_pos scale < 1, 1, scaleZ > } 
    }
    // Negative part:
    object{
      scene_mesh scale < 1, 1, scaleZ >
      clipped_by{ clip_neg }
      texture { tx_neg scale < 1, 1, scaleZ > } 
    }
    // Null part:
    object{ 
      scene_mesh scale < 1, 1, scaleZ >
      clipped_by{ clip_zer } 
      texture { tx_zer } 
    }
  #else
    object{ scene_mesh translate <-(maxX)/2.0,-(maxY)/2.0,0 > scale < 1, 1, scaleZ > } 
//     sphere{ <0,0,0>, 0.5*radiusXY }
    texture { tx_pos scale < 1, 1, scaleZ > }
  #end
}

union{
  light_source { <20000, 5000, 1000> color rgb 0.80*<1,1,1> }
  light_source { <20000, -5000, -1000> color rgb 0.30*<1,1,1> }
  rotate -(30 + viewElev/2)*y
  rotate -viewAzim*z
}

#declare camera_dir = vrotate(vrotate(<1,0,0>, -viewElev*y), -viewAzim*z);
#declare camera_disp = viewDist*camera_dir;
#declare camera_radius = 0.95*radiusXY;
#declare camera_focus = <0,0,0>;
#declare image_rel_width = image_width;
#declare image_rel_height = image_height;

#include "camera.inc"
