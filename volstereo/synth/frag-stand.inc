// Last edited on 2003-04-17 09:03:10 by stolfi

// Grid dimensions: 
#declare grid_size = < 180.0, 230, 0 >;

// Half-width of grid lines:
#declare grid_line_radius = 0.2;

// Paper dimensions:
#declare frame_size = < 210.0, 297.0, 0.1 >;

// Corners of paper frame:
#declare frame_min = < -19.0, -27.0, - frame_size.z >; 
#declare frame_max = frame_min + frame_size; 

// Gap between frame and mirror:
#declare frame_gap = 0.1;

// Mirror dimensions (similar to US Letter):
#declare mirror_size = < 210.0, 297.0, 2.5 >;

// Corners of mirror:
#declare mirror_min = < -19.0, -27.0, frame_min.z - frame_gap - mirror_size.z >; 
#declare mirror_max = mirror_min + mirror_size; 

// Table top:
#declare table_Z = mirror_min.z - eps;


#declare tx_paper = 
  texture {
    pigment { color rgb < 1.000, 1.000, 1.000 > }
    finish { ambient 0.0 diffuse 1.0 }
  }
  
#declare tx_black_ink =
  texture {
    pigment { color rgb < 0.020, 0.020, 0.020 > }
    finish { ambient 0.0 diffuse 0.9 specular 0.1 roughness 0.005 }
  }
  
#declare tx_glass =
  texture {
    pigment { color rgb < 0.980, 0.980, 0.980 > filter 0.980 }
    finish { 
      ambient 0.0 
      diffuse 0.010 
      specular 0.005 roughness 0.002 
      reflection 0.050
    }
  }
      
#declare in_glass =
  interior {
    ior 1.5
  }

#declare tx_mirror =
  texture {
    pigment { color rgb < 1.000, 1.000, 1.000 > }
    finish { 
      ambient 0.0 
      diffuse 0.005 
      reflection 0.980
    }
  }

#declare stand_mirror_glass = 
  box { 
    mirror_min - eps * z,
    mirror_min + mirror_size
    texture { tx_glass }
    interior { in_glass }
  }

#declare stand_mirror_metal = 
  plane { 
    -z, 0
    texture { tx_mirror }
    translate mirror_min.z * z
  }

#declare stand_mirror =
  intersection { 
    object { stand_mirror_glass }
    object { stand_mirror_metal }
    no_shadow
  }

#declare stand_frame_paper =
  box { 
    frame_min,
    frame_max
  }
  
#declare stand_frame_hole = 
  difference { 
    box { 
      < +005.0, +005.0, -010.0 > 
      < grid_size.x - 005.0, grid_size.y - 005.0, +010.0 >
    }
    box { 
      < +060.0, 0000.0, -011.0 > 
      < +120.0, +020.0, +011.0 > 
    }
    box { 
      < 0000.0, +170.0, -011.0 > 
      < +020.0, +230.0, +011.0 > 
    }
    box { 
      < +160.0, +170.0, -011.0 > 
      < +180.0, +230.0, +011.0 > 
    }
    box { 
      < +040.0, +220.0, -011.0 > 
      < +060.0, +230.0, +011.0 > 
    }
    box { 
      < +120.0, +220.0, -011.0 > 
      < +140.0, +230.0, +011.0 > 
    }
  }
  
#declare stand_frame_grid = 
  union {
    #declare dx = 0;
    #while (dx <= grid_size.x)
      cylinder { < dx, 0, 0 > < dx, grid_size.y, 0 >, grid_line_radius }
      #declare dx = dx + 10.0;
    #end
    #declare dy = 0;
    #while (dy <= grid_size.y)
      cylinder { < 0, dy, 0 > < grid_size.x, dy, 0 >, grid_line_radius }
      #declare dy = dy + 10.0;
    #end
    sphere { < 0, 0, 0 >, grid_line_radius }
    sphere { < grid_size.x, 0, 0 >, grid_line_radius }
    sphere { < 0, grid_size.y, 0 >, grid_line_radius }
    sphere { < grid_size.x, grid_size.y, 0 >, grid_line_radius }
    scale < 1, 1, 0.1 >
    texture { tx_black_ink }
  }

#declare stand_frame =
  union{
    difference { 
      object { stand_frame_paper }
      object { stand_frame_hole }
      texture { tx_paper }
    }
    difference { 
      object { stand_frame_grid translate frame_size.z * z }
      object { stand_frame_hole }
      texture { tx_black_ink }
    }
  }
  
#declare stand = 
  union {
    object { stand_mirror } 
    object { stand_frame }
  }

