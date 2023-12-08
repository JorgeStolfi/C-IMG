// Last edited on 2003-04-17 12:44:56 by stolfi

#declare light_color = <0.848, 0.792, 0.775>;
#declare light_pwr = 0.50;

#declare light_elev1 = -50;
#declare light_azim1 = 135;

#declare light_elev2 = -30;
#declare light_azim2 = -167;

#declare light_elev3 = -30;
#declare light_azim3 = -13;

#declare light_dir1 = vrotate(vrotate(<1,0,0>, light_elev1 * y), light_azim1*z);
#declare light_dir2 = vrotate(vrotate(<1,0,0>, light_elev2 * y), light_azim2*z);
#declare light_dir3 = vrotate(vrotate(<1,0,0>, light_elev3 * y), light_azim3*z);

#declare light_ctr = < +108.0, +140.0, 0000.0 >;
#declare light_dist = 600.0;

#declare lights = 
  union {
    light_source {
      light_ctr + light_dist * light_dir1
      color rgb light_pwr * light_color
    }

    light_source {
      light_ctr + light_dist * light_dir2
      color rgb light_pwr * light_color
    }

    light_source {
      light_ctr + light_dist * light_dir3
      color rgb light_pwr * light_color
    }
  }
  
// #declare tx_light_vec = 
//   texture { pigment { color rgb <0,1,1> } finish { ambient 0.5 diffuse 0.5 } } 
//   
// cylinder { light_ctr, light_ctr + 60 * light_dir1, 2.0 texture { tx_light_vec } } 
// cylinder { light_ctr, light_ctr + 60 * light_dir2, 2.0 texture { tx_light_vec } } 
// cylinder { light_ctr, light_ctr + 60 * light_dir3, 2.0 texture { tx_light_vec } } 
