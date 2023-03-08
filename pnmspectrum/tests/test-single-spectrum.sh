#! /bin/bash
# Last edited on 2010-10-04 12:04:05 by stolfi

name="$1"; shift

infile="data/${name}.ppm"
outppm="out/spectrum-${name}.ppm"
outtbl="out/spectrum-${name}.txt"

../pnmspectrum \
     -maxPower 4.0 \
     -outputImage ${outppm} \
       -center \
       -scale log 1.0e-3 \
       -zeroMean \
     -outputTable ${outtbl} -exact \
     -verbose \
     ${infile}
display ${outppm} ${infile}
