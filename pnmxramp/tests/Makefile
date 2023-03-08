# Last edited on 2012-12-21 18:13:52 by stolfilocal

PROG := pnmxramp
PROGDIR := ..

TESTMAPS := \
  out/RC-0640-0480.ppm \
  out/RGBC-0640-0480.ppm \
  out/GRBM-0640-0480.ppm \
  out/BRGY-0640-0480.ppm \
  out/BRGY-0640-0480.ppm \
   \
  out/VL-0640-0480.pgm \
  out/HL-0640-0480.pgm \
  out/DL-0640-0480.pgm \

all: ${TESTMAPS}

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

out/RC-0640-0480.ppm: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -size 0640 0480 \
          -channels 3 \
          -maxval 100 \
          0054 0112  0.000 0.000 0.000 \
          0566 0112  1.000 0.000 0.000 \
          0054 0368  0.000 1.000 1.000 \
          > $@
	display -gamma 2.200 $@
          
out/RGBC-0640-0480.ppm: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -size 0640 0480 \
          0054 0112  1.000 0.000 0.000 \
          0566 0112  0.000 1.000 0.000 \
          0054 0368  0.000 0.000 1.000 \
          > $@
	display -gamma 2.200 $@
        
out/GRBM-0640-0480.ppm: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -size 0640 0480 \
          0054 0112  0.000 1.000 0.000 \
          0566 0112  1.000 0.000 0.000 \
          0054 0368  0.000 0.000 1.000 \
          > $@
	display -gamma 2.200 $@

out/BRGY-0640-0480.ppm: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -size 0640 0480 \
          0054 0112  0.000 0.000 1.000 \
          0566 0112  1.000 0.000 0.000 \
          0054 0368  0.000 1.000 0.000 \
          > $@
	display -gamma 2.200 $@

out/VL-0640-0480.pgm: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -size 0640 0480 \
          0054 0112  0.900 0.900 0.900 \
          0566 0112  0.900 0.900 0.900 \
          0054 0368  0.100 0.100 0.100 \
          > $@
	display -gamma 2.200 $@

out/HL-0640-0480.pgm: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -size 0640 0480 \
          -channels 1\
          0054 0112  0.900 0.900 0.900 \
          0566 0112  0.100 0.100 0.100 \
          0054 0368  0.900 0.900 0.900 \
          > $@
	display -gamma 2.200 $@

out/DL-0640-0480.pgm: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
          -size 0640 0480 \
          -channels 1\
          0054 0112  0.100 0.100 0.100 \
          0566 0112  0.900 0.900 0.900 \
          0054 0368  0.900 0.900 0.900 \
          > $@
	display -gamma 2.200 $@

clean::
	rm -f ${TESTMAPS}

