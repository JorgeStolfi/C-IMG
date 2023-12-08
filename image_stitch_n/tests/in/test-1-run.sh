#! /bin/bash
# Last edited on 2012-01-04 00:35:16 by stolfi

${PROGDIR}/${PROG} \
  -image A 0 0  120 0  120 240  0 240 \
  -image B 0 0  120 0  120 240  0 240 \
  -image C 0 0  120 0  120 240  0 240 \
  -image D 0 0  120 0  120 240  0 240 \
  \
  -match A B in/test-1-A-B-pairs.txt \
  -match A C in/test-1-A-C-pairs.txt \
  -match C D in/test-1-C-D-pairs.txt \
  -match B D in/test-1-B-D-pairs.txt \
  -maxIter 5 \
  -maxErr 0 \
  -outPrefix out/test-1 \
  -verbose YES
  
