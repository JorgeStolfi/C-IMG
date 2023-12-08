# Last edited on 2009-01-06 17:51:48 by stolfi

PROG := pnmpairtopng
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  butt-w
  
NOMTESTS := \
  ramp \
  ic-logo
  
all:
	for t in ${TESTS}; do \
          ${MAKE} TEST=$$t single ; \
        done
        
clean::
	for t in ${TESTS} ${NONTESTS}; do \
          ${MAKE} TEST=$$t clean-single ; \
        done

SHOW := display -title '%f'

######################################################################
# Sub-make - caller must define ${TEST}

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

OLDPROG := ${STOLFIHOME}/bin/linear-ppm-pair-to-png

single: ${TEST}.png ${TEST}-a.pgm ${TEST}-old-a.pgm ${TEST}-a-dif.pgm
	${SHOW} -filter Point -geometry '200%' \
          ${TEST}-0.ppm ${TEST}-1.ppm ${TEST}.png \
          ${TEST}-a.pgm ${TEST}-old-a.pgm \
          ${TEST}-a-dif.pgm

${TEST}.png: ${PROGDIR}/${PROG} ${TEST}-0.ppm ${TEST}-1.ppm
	${PROGDIR}/${PROG} ${TEST}-0.ppm ${TEST}-1.ppm > ${TEST}.png
	identify ${TEST}.png

${TEST}-old.png: ${TEST}-0.ppm ${TEST}-1.ppm ${OLDPROG}
	${OLDPROG} ${TEST}-0.ppm ${TEST}-1.ppm > ${TEST}-old.png
	identify ${TEST}.png

%-a.pgm: %.png
	convert $*.png -channel Matte $*-a.pgm

${TEST}-a-dif.pgm: ${TEST}-a.pgm ${TEST}-old-a.pgm
	pnmxarith \
            -scale 64 -offset 127 \
            -subtract ${TEST}-a.pgm ${TEST}-old-a.pgm \
          > ${TEST}-a-dif.pgm

clean-single:
	-rm -f ${TEST}.png ${TEST}-old.png 
	-rm -f ${TEST}-a.pgm ${TEST}-old-a.pgm ${TEST}-a-dif.pgm

endif
# end ${TEST} sub-make.
######################################################################

######################################################################
# For the "ramp" test:

MAKERAMP := ppmxramp

ramp-0.ppm: Makefile ${MAKERAMP}
	${MAKERAMP} \
            256 256 \
            000 000 rgbi:0/0/0 \
            255 000 rgbi:1/1/1 \
            255 255 rgbi:1/1/1 \
          > ramp-0.ppm
	${SHOW} ramp-0.ppm

ramp-1.ppm: Makefile ${MAKERAMP}
	${MAKERAMP} \
            256 256 \
            000 000 rgbi:0/0/0 \
            000 255 rgbi:1/1/1 \
            255 255 rgbi:1/1/1 \
          > ramp-1.ppm
	${SHOW} ramp-1.ppm
          
# end creation of images for "ramp" test.
######################################################################
