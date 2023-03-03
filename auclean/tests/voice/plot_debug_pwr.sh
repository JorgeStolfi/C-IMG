#! /bin/bash
# Last edited on 2023-03-02 18:00:58 by stolfi

dfile="$1"; shift  # File name with frame debug data.
fmax="$1"; shift   # Max freq to plot.

tmp="/tmp/$$"

tfile="${tmp}.png"
pfile="${dfile/.pwr/.png}"

export GDFONTPATH=${STOLFILOCAL}/ttf

gnuplot << EOF
set term png size 3000,1900 font "arial,14"
set output "${tfile}"

set multiplot layout 5,1 title "Noise removal iteration"
set nokey
set xrange [-10:${fmax}]
# ----------------------------------------------------------------------
# Input frame power
set logscale y
set yrange [0.00001:]
set format y "%12.8f"
set ylabel "input"
plot "${dfile}" using 1:2 with linespoints pt 7 ps 1.0 lc rgb '#008800'
# ----------------------------------------------------------------------
# Assumed input noise power
set logscale y
set yrange [0.00001:]
set format y "%12.8f"
set ylabel "noise"
plot "${dfile}" using 1:3 with linespoints pt 7 ps 1.0 lc rgb '#0000aa'
# ----------------------------------------------------------------------
# Frequency weights for least squares fitting
set nologscale y
set yrange [0.001:]
set format y "%12.4f"
set ylabel "weight"
plot "${dfile}" using 1:4 with linespoints pt 7 ps 1.0 lc rgb '#aa0000'
# ----------------------------------------------------------------------
# Fitted noise power
set logscale y
set yrange [0.00001:]
set format y "%12.8f"
set ylabel "fit noise"
plot "${dfile}" using 1:5 with linespoints pt 7 ps 1.0 lc rgb '#999999
# ----------------------------------------------------------------------
# Audio power
set logscale y
set yrange [0.00001:]
set format y "%12.8f"
set ylabel "audio"
plot "${dfile}" using 1:6 with linespoints pt 7 ps 1.0 lc rgb '#008800'
# ----------------------------------------------------------------------
unset multiplot

pause mouse

EOF

if [[ -s ${tfile} ]]; then
  convert ${tfile} -resize '50%' ${pfile}
  display -loop 0 -delay 100x1 ${pfile}
  rm ${tfile}
fi

