#! /bin/bash
# Last edited on 2024-10-17 12:32:30 by stolfi

src="../../../../../JSLIBS/libimg/tests/330_test_images/out"
dst="./"

for tags in 1024x1024-1:.pgm 0320x0240-3:-sma.ppm ; do 
  tag_src="${tags/:*/}"
  tag_dst="${tags/*:/}"
  for img in chopsea grittie bullsex bullsqr ; do
    convert ${src}/test-${tag_src}-${img}.png -flip ${img}${tag_dst}
  done
done
