#! /bin/bash
# Last edited on 2017-07-28 02:48:46 by jstolfi

cmd="${0##*/}"
usage="${cmd} IMAGENAME"

# Generates PostScript files containing plots of the color calibration
# data which is read from files "in/{IMAGENAME}-med.txt" and
# "out/{IMAGENAME}-plo.txt"

if [ $# -ne 1 ]; then
  echo "usage: ${usage}"; exit 1
fi

img="$1"; shift;
echo "plotting in/${img}-med.txt -> out/${img}-med.ps"

gnuplot <<EOF
set terminal postscript
set output "out/${img}-med.ps"
set xlabel "yi*"
set xrange [0:1]
set ylabel "yi',yi''"
set yrange [0:255]
# Plots RGB' and RGB'' as a function of RGB*. 
plot \
  "in/${img}-med.txt" using (0.1**column(2)):3 title "R'" w l lt 3 lw 3, \
  "in/${img}-med.txt" using (0.1**column(2)):4 title "G'" w l lt 1 lw 3, \
  "in/${img}-med.txt" using (0.1**column(2)):5 title "B'" w l lt 2 lw 3, \
  \
  "in/${img}-med.txt" using (0.1**column(2)):6 title "R''" w l lt 3 lw 3, \
  "in/${img}-med.txt" using (0.1**column(2)):7 title "G''" w l lt 1 lw 3, \
  "in/${img}-med.txt" using (0.1**column(2)):8 title "B''" w l lt 2 lw 3
quit
EOF

atril out/${img}-med.ps

echo "plotting out/${img}-plo.txt -> out/${img}-plo.ps"

gnuplot <<EOF
set terminal postscript
set output "${img}-plo.ps"
set xlabel "log(yi*)"
set xrange [-2:0]
set ylabel "log(yi'/yi'')"
set yrange [-2.8:0.2]
plot \
  "${img}-plo.txt" using 1:4 title "R" w l lt 3 lw 3, \
  "${img}-plo.txt" using 2:5 title "G" w l lt 1 lw 3, \
  "${img}-plo.txt" using 3:6 title "B" w l lt 2 lw 3
quit
EOF

atril ${img}-plo.ps
