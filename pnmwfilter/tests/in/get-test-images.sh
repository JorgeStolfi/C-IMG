#! /bin/bash
# Last edited on 2010-08-05 00:22:02 by stolfi

ib2="${STOLFIHOME}/projects/imgbank/imgbank2"

# convert ${ib2}/000/000/067/010.png -resize '075%' -crop '256x256+060+000' i0010.pgm
# convert ${ib2}/000/000/104/037.png -resize '100%' -crop '256x256+128+128' i0037.pgm
# convert ${ib2}/000/000/039/045.png -resize '100%' -crop '256x256+080+030' i0045.pgm
# convert ${ib2}/000/001/006/000.png -resize '100%' -crop '256x256+128+128' i1000.pgm
# convert ${ib2}/000/001/002/002.png -resize '100%' -crop '256x256+128+128' i1002.pgm
# convert ${ib2}/000/001/018/006.png -resize '075%' -crop '256x256+080+000' i1006.pgm
# convert ${ib2}/000/001/018/009.png -resize '100%' -crop '256x256+012+060' i1009.pgm

vms="${STOLFIHOME}/projects/voynich/docs/RZCD/cdrom-copy"

# convert ${vms}/fol116v.bmp -resize '100%' -crop '256x256+166+102' i2001.pgm

prg="${STOLFIHOME}/programs/c"

convert ${prg}/COURSES/mc919-fi/IN/fruta.ppm                     -resize '100%' -crop '256x256+000+000' i3001.ppm
convert ${prg}/IMG/pnmfield/tests/pnmadjust/data/jaburu-1-o.jpg  -resize '100%' -crop '256x256+020+000' i3002.ppm
convert ${prg}/IMG/pnmift/tests/data/footb.ppm                   -resize '100%' -crop '256x256+000+000' i3003.ppm
convert ${prg}/IMG/volstereo/tests/m21.ppm                       -resize '100%' -crop '256x256+030+300' i3004.ppm
convert ${ib2}/000/000/039/726.png                               -resize '100%' -crop '256x256+000+240' i3005.ppm
convert ${ib2}/000/000/040/144.png                               -resize '100%' -crop '256x256+080+160' i3006.ppm
convert ${ib2}/000/000/070/004.png                               -resize '100%' -crop '256x256+010+240' i3007.ppm
convert ${ib2}/000/001/000/012.png                               -resize '100%' -crop '256x256+020+100' i3008.ppm
convert ${ib2}/000/001/062/007.png                               -resize '100%' -crop '256x256+000+000' i3009.ppm

convert i0???.pgm +append strip0.png
convert i1???.pgm +append strip1.png
convert i2???.pgm +append strip2.png
convert i3???.ppm +append strip3.png
convert strip?.png -append strip.png
display strip.png
rm -f strip?.png

