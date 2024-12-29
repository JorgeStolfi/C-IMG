# Last edited on 2008-01-14 19:35:52 by stolfi

PROG := geostereo
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: eye-disp.pgm

eye-disp.pgm: eye-0.pgm eye-1.pgm ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -window 3 3 -nscales 0 \
	    -mindisp -20 -maxdisp +20 \
	    -output eye \
	    eye-0.pgm eye-1.pgm
	display eye-0.pgm & display eye-1.pgm & display eye-disp.pgm eye-score.pgm

eye-%.pgm: eye-%.ppm
	ppmtopgm < eye-$*.ppm > eye-$*.pgm

eye-%.ppm: eye-%.pov scene.inc camera.inc do-run-povray
	do-run-povray $* -preview -good
	cp -v `cat .last`/eye-$*.ppm eye-$*.ppm
