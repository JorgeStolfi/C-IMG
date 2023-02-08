#! /bin/bash
# Last edited on 2023-02-08 09:33:22 by stolfi

ppmfile=$1; shift;
hstfile="${ppmfile:r}.hst"
gbhfile="${ppmfile:r}.gbh"

pnmsmooth -size 3 3 ${ppmfile} \
  | ppmhist - \
  > ${hstfile}

cat ${hstfile} \
  | gawk \
      -v rs=032 -v gs=012 -v bs=007 \
      -v r0=255 -v g0=000 -v b0=000 \
      ' /^[ ]*[0-9]/{ \
          r = $1-rs; g = $2-gs; b = $3-bs; n = $5; \
          w = r; if (w+0 == 0) { w = 1; } \
          printf "%6.3f %6.3f %6.3f %7d\n", \
            (r-r0)/w, g/w*r0 - g0, b/w*r0 - b0, n; \
        } \
      ' \
  > ${gbhfile}

#   -v rs=000 -v gs=000 -v bs=000 \
#   

gnuplot <<EOF
# set xrange [-1.5:+1.5]
# set yrange [-1.5:+1.5]
set xrange [050:300]
set yrange [050:300]
set size ratio -1 1.0000,1.0000
plot "${gbhfile}" using 2:3 with points
pause 300
EOF
