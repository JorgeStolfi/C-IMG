#! /usr/bin/gawk -f
# Last edited on 2010-06-17 23:35:18 by stolfi

BEGIN{ 
  
  for (y = 0; y < 1000; y++)
    { for (x = 0; x < 1000; x++) 
        { bx = (x % 2);
          by = (y % 2);
          if (x < 500) { m = 20000; } else { m = 15000; }
          val = (2*by + bx)*m + (x + y);
          printf "%c%c", val/256, val%256; 
        }
    }
    
  for (z = 0; z < 1472; z++)
    { val = z;
      printf "%c%c", val/256, val%256; 
    }
}

