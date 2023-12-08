#! /bin/csh -f
# Last edited on 2000-06-23 04:06:52 by stolfi

# set infile = small.ppm
# set col = 4
# set row = 0
# set wd = 1
# set ht = 9
# set shift = ( 00.000 -0.500  00.000 00.000  00.000 +0.333 )

set infile = test.ppm
set col = 66
set row = 140
set wd = 1
set ht = 50
set shift = ( 0.000 -0.500  0.000 0.000  0.000 0.000 )

plot-image-slice ${col} ${row} ${wd} ${ht} < ${infile} > ${infile:r}-slice.gif
../pnmbend -shift ${shift} ${infile} > res.ppm
if (! ( -z res.ppm ) ) then
  plot-image-slice ${col} ${row} ${wd} ${ht} < res.ppm > res-slice.gif
  xv res.ppm ${infile} & 
endif
