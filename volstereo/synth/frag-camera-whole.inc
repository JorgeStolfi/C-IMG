// Last edited on 2003-04-17 09:03:59 by stolfi

#declare camera_ctr = < +093.0, +123.0, 0000.0 >;
#declare camera_dir = vrotate(-y, -32*x);
#declare camera_dist = 605.0;
#declare camera_size = 0.720;

camera { 
  location camera_ctr + camera_dist * camera_dir
  right -0.512 * camera_size * x 
  up    +0.384 * camera_size * y
  sky   z
  look_at camera_ctr
}
  
// #declare tx_camera_vec = 
//   texture { pigment { color rgb <1,0.8,1> } finish { ambient 0.5 diffuse 0.5 } } 
//   
// #declare camera_p1 = camera_ctr + 60 * camera_dir;
// #declare camera_p2 = camera_p1 + 7*camera_dir;
// cylinder { camera_ctr, camera_p1, 2.0 texture { tx_camera_vec } } 
// cone { camera_p1, 4.0, camera_p2, 0.0 texture { tx_camera_vec } } 
