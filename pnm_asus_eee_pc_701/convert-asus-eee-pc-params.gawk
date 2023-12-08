#! /usr/bin/gawk -f
# Last edited on 2008-11-24 00:12:51 by stolfi

BEGIN {
  WHITE = 255.0;
  BLACK = 16.0;
  printf "  \"   {bbb} {BRGHT} \\n\" \\\n";
  printf "  \"   ----- ------- \\n\" \\\n";
  for (bbb = 25; bbb <= 250; bbb += 25)
    { b = (bbb - 128 - BLACK)/(WHITE-BLACK);
      printf "  \"    %03d  %+7.4f \\n\" \\\n", bbb, b;
    }
  printf "  \"\\n\" \\\n";
  printf "  \"   {ccc} {CTRST} \\n\" \\\n";
  printf "  \"   ----- ------- \\n\" \\\n";
  for (ccc = 6; ccc <= 30; ccc += 6)
    { c = ccc/16.0/(WHITE-BLACK);
      printf "  \"    %03d  %+7.4f \\n\" \\\n", ccc, c;
    }
}
