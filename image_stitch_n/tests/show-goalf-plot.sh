#! /bin/bash
# Last edited on 2011-12-30 21:36:23 by stolfilocal

gnuplot <<EOF
set term x11
plot \
  "out/plot.txt" using 2:3  with linespoints pt 7 lt 1 lc rgb '#ff0000', \
  "out/plot.txt" using 2:4  with linespoints pt 7 lt 1 lc rgb '#aa5500', \
  "out/plot.txt" using 2:5  with linespoints pt 7 lt 1 lc rgb '#226600', \
  "out/plot.txt" using 2:6  with linespoints pt 7 lt 1 lc rgb '#00aa00', \
  "out/plot.txt" using 2:7  with linespoints pt 7 lt 1 lc rgb '#008888', \
  "out/plot.txt" using 2:8  with linespoints pt 7 lt 1 lc rgb '#0033ff', \
  "out/plot.txt" using 2:9  with linespoints pt 7 lt 1 lc rgb '#7700ff', \
  "out/plot.txt" using 2:10 with linespoints pt 7 lt 1 lc rgb '#aa0044'
pause 200
quit
EOF
