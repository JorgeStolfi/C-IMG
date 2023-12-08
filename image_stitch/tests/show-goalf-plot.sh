#! /bin/bash
# Last edited on 2012-01-23 17:57:49 by stolfilocal

testnum="$1"; shift
infile="$1"; shift

gnuplot <<EOF
set term x11
set title "Test ${testnum} - Goal function"
plot \
  "${infile}" using 2:3  with linespoints pt 7 lt 1 lc rgb '#ff0000', \
  "${infile}" using 2:4  with linespoints pt 7 lt 1 lc rgb '#aa5500', \
  "${infile}" using 2:5  with linespoints pt 7 lt 1 lc rgb '#226600', \
  "${infile}" using 2:6  with linespoints pt 7 lt 1 lc rgb '#00aa00', \
  "${infile}" using 2:7  with linespoints pt 7 lt 1 lc rgb '#008888', \
  "${infile}" using 2:8  with linespoints pt 7 lt 1 lc rgb '#0033ff', \
  "${infile}" using 2:9  with linespoints pt 7 lt 1 lc rgb '#7700ff', \
  "${infile}" using 2:10 with linespoints pt 7 lt 1 lc rgb '#aa0044'
pause 200
quit
EOF
