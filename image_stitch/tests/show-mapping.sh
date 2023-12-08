#! /bin/bash
# Last edited on 2012-01-23 18:48:27 by stolfilocal

# Arguments are the test number, the mapped outline files of the two images,
# the mapped point pairs file, and a magnification factor.

testnum="$1"; shift
dom1file="$1"; shift
dom2file="$1"; shift
pairfile="$1"; shift
mag="$1"; shift;

tmp="/tmp/$$"
tmpimage="${tmp}-t.png"
outimage="${pairfile%%.*}.png"
  
export GDFONTPATH="./ttf"

gnuplot <<EOF
set term png truecolor size 1200,1200 font "arial,18"
set output "${tmpimage}"
# set term x11
set size ratio -1
set title "Test ${testnum} - mapped outlines and points"
unset key
mag=${mag}
mdx(dummy) = (column(3)+column(7))/2
mdy(dummy) = (column(4)+column(8))/2
dpx(dummy) = mag*(column(7)-column(3))
dpy(dummy) = mag*(column(8)-column(4))

set style arrow 1  lw 3 lc rgb '#00ccff' head filled size character 1,20,90

plot \
  "${dom1file}" using 1:2 with lines lt 1 lw 3 lc rgb '#ff0088', \
  "${dom2file}" using 1:2 with lines lt 1 lw 3 lc rgb '#00aa00', \
  "${pairfile}" using (mdx(0)):(mdy(0)):(dpx(0)):(dpy(0)) with vectors arrowstyle 1, \
  "" using (mdx(0)):(mdy(0)) with points lt 0 pt 7 ps 2.0
# pause 200
quit
EOF

convert ${tmpimage} -resize '50%' ${outimage}
display ${outimage}
rm -fv ${tmpimage}
