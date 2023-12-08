#! /usr/bin/gawk -f
# Last edited on 2000-06-23 03:53:37 by stolfi

BEGIN{
  r = 4;
  wd = 2*r+1; ht = 2*r+1;
  printf "P3\n";
  printf "%d %d\n", wd, ht;
  printf "255\n";
  for(y=-4;y<=4;y++){
    for(x=-4;x<=4;x++){
      r = f(x, y+0.500); g = f(x,y); b = f(x, y-0.333);
      printf "%03d %03d %03d  ", r,g,b;
    }
    printf "\n";
  }
}

function f(x,y){
  return int(50+150*exp(-(x^2+y^2)/2) + 0.5);
}
