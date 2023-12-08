#! /bin/bash
# Last edited on 2018-09-10 18:42:38 by stolfilocal

# Usage:
#
#   plot_data_points.sh {SHOW} {DIR}
#
# where
#
#   {SHOW} "SHOW" to display, "NOSHOW" not to.
#   {DIR} is the file name prefix.
#
# Reads {DIR}/data.txt, which is supposed to contain the data points for 
# linear regression, one per line, with blank lines between 
# data points of different window positions. Writes {DIR}/data.png, the
# plot file.

show="$1"; shift
dir="$1"; shift

tmp=/tmp/$$

dfile="${dir}/data.txt"
tfile="${tmp}.txt"
pfile="${dir}/data.png"

# Create a random sample of the input data:
cat ${dfile} \
  | gawk \
      ' /[0-9]/ { 
          p = $1;
          d = $2;
          f = $5;
          e2 = $7;
          m2 = $8;
          s2 = $9;
          un = $10;

          if ((d % 100) != 17) { next; }
          if (f == 0) { printf "\n"; }
          printf "%5d %5d %6.2f %18.8f\n", d, f, un, s2;
        }
      ' \
  > ${tfile}

export GDFONTPATH=.

gnuplot <<EOF
set term png size 2800,1000 font "arial,20"
set xtics out 5; set mxtics 5
set grid xtics lt 1 lw 3 lc rgb '#ffddaa', lt 1 lw 1.5 lc rgb '#ffddaa'
set grid mxtics
set output "${tmp}.png"
set yrange [-0.02:+1.02]
unset key
plot "${tfile}" using 3:4 title 's2f' with linespoints pt 7 lc rgb '#008800'
quit
EOF

convert ${tmp}.png -resize '50%' ${pfile}

if [[ ( "/${show}" == "/SHOW" ) && ( -s ${pfile} ) ]]; then
  display -title '%f' ${pfile}
fi

rm -fv ${tmp}.*
