#! /bin/bash
# Last edited on 2024-12-25 22:17:35 by stolfi

# Plots a test histogram produced by "pnmxhist"

hisfile="$1"; shift
pngfile="$1"; shift
maxval="$1"; shift

tmp="/tmp/$$"

gnuplot <<EOF
  set term png truecolor rounded medium size 400,400
  set output "${tmp}-count-abs.png"
  set title "count (abs)"
  set yrange [-1: ]
  set xrange [-1:(${maxval}+1)]
  set nokey
  plot "${hisfile}" using 1:2 with histeps
  quit
EOF

gnuplot <<EOF
  set term png truecolor rounded medium size 400,400
  set output "${tmp}-count-rel.png"
  set title "count (rel)"
  set yrange [-0.000001: ]
  set xrange [-1:(${maxval}+1)]
  set nokey
  plot "${hisfile}" using 1:3 with histeps
  quit
EOF

gnuplot <<EOF
  set term png truecolor rounded medium size 400,400
  set output "${tmp}-accum-abs.png"
  set title "accum (abs)"
  set yrange [-1: ]
  set xrange [-1:(${maxval}+1)]
  set nokey
  plot \
    "${hisfile}" using 1:4 with histeps, \
    "${hisfile}" using 1:6 with histeps
  quit
EOF

gnuplot <<EOF
  set term png truecolor rounded medium size 400,400
  set output "${tmp}-accum-rel.png"
  set title "accum (rel)"
  set yrange [-0.000001:1.000001]
  set xrange [-1:(${maxval}+1)]
  set nokey
  plot \
    "${hisfile}" using 1:5 with histeps, \
    "${hisfile}" using 1:7 with histeps
  quit
EOF

convert \
  +append \
  ${tmp}-count-abs.png \
  ${tmp}-count-rel.png \
  ${tmp}-accum-abs.png \
  ${tmp}-accum-rel.png \
  ${pngfile}

rm -f ${tmp}-*.png
