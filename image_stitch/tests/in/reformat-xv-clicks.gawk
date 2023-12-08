#! /usr/bin/gawk -f 

($3 == "=") {  
  x = $1; gsub(/[,]$/, "", x);
  y = $2;
  printf "( %5d %5d )\n", x, y; next;
}

//{ print; next; }
