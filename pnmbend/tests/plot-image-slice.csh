#! /bin/csh -f
# Last edited on 2008-02-04 20:07:49 by stolfi

set usage = "$0 COL ROW WD HT < INFILE.ppm > PLOT.gif"

# Plots R, G, B along pixels of an image

if ($#argv != 4) then
  echo "usage: ${usage}"; exit 1
endif

set row = "$1"; shift;
set col = "$1"; shift;
set wd = "$1"; shift;
set ht = "$1"; shift;

set tmp = "/tmp/$$"
set infile = "${tmp}.pixs"
set bmfile = "${tmp}.ppm"
set gffile = "${tmp}.gif"

pnmcut ${row} ${col} ${wd} ${ht} \
  | pnmnoraw \
  | sed -e '/^[#]/d' \
  | tr ' ' '\012' \
  | gawk \
      ' BEGIN{n=-4;} \
        /./{ \
          s[n+0]=$1; n++; \
          if(n==3){printf "%3d %3d %3d\n",s[0],s[1],s[2];n=0;} \
        } ' \
  | gawk \
      -v diff=1 \
      ' ((FNR==1)&&(diff==1)){r=$1; g=$2; b=$3; next;} \
        /./ { \
          printf "%4d %4d %4d\n", $1-r, $2-g, $3-b; \
          if(diff==1){r=$1; g=$2; b=$3}; \
        } ' \
  >  ${infile}

gnuplot <<EOF
set terminal pbm color medium
set output "${bmfile}"
set xlabel "pixel"
set ylabel "intensity"
# set yrange [-5:+261]
# set yrange [-261:+261]
set size 1.5,2.0
set xzeroaxis
plot \
  "${infile}" using 1 title "red" with linespoints lt 1 pt 1, \
  "${infile}" using 2 title "grn" with linespoints lt 2 pt 1, \
  "${infile}" using 3 title "blu" with linespoints lt 3 pt 1
quit
EOF

if ( ( -r ${bmfile} ) && ( ! ( -z ${bmfile} ) ) ) then
  ppmtogif < ${bmfile} > ${gffile}
  cat ${gffile}
  ( xv ${gffile}; /bin/rm -f  ${gffile} ) &
endif

/bin/rm -f ${bmfile}
 
