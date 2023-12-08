#! /bin/bash
# Last edited on 2012-01-04 00:35:04 by stolfi

${PROGDIR}/${PROG} \
  -image A 0 0  200 0  200 200  0 200 \
  -image B 0 0  200 0  200 200  0 200 \
  -image C 0 0  200 0  200 200  0 200 \
  -image D 0 0  200 0  200 200  0 200 \
  \
  -match A B in/test-2-A-B-pairs.txt \
  -match A C in/test-2-A-C-pairs.txt \
  -match C D in/test-2-C-D-pairs.txt \
  -match B D in/test-2-B-D-pairs.txt \
  -maxIter 5 \
  -maxErr 0 \
  -outPrefix out/test-2 \
  -eqRotation YES \
  -center 150 150 \
  -verbose YES
  
