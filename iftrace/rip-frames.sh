#! /bin/bash 
# Last edited on 2010-06-09 01:38:01 by stolfilocal

vid=$1; shift

mkdir -p ${vid}/frames

for i in 0 10 20 30 40 50 60 70 80 90 100 110 120 130 ; do
  ix=`printf '%05d' ${i}`
  mkdir -p ${vid}/frames/${ix}
  convert ${vid}/video.mpg PPM:- \
    | convert "PPM:-[${i}]"  ${vid}/frames/${ix}/frame.ppm
  ls -l ${vid}/frames/${ix}/frame.ppm
done

( cd ${vid}/frames && display -title '%d' */frame.ppm )
