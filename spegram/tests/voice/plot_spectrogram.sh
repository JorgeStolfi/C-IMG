#! /bin/bash
# Last edited on 2006-10-27 08:26:16 by stolfi

datafile="$1"

gnuplot <<EOF
set terminal X11
set hidden3d
splot \
  "< egrep -v '^[ ]*[\#]' ${datafile}" using 1:2:3 with lines
pause 300
EOF
