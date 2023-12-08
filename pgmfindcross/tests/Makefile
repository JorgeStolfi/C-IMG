# Last edited on 2009-01-06 17:49:08 by stolfi

PROG := pgmfindcross
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

.SUFFIXES:

SHOW := display -filter Point

TESTS := \
  f0ay:-darklines@-linewd@2:5:0.33333

NOTTESTS := \
  fxx:-darklines:5:1.0000 \
  fyy:@:5:1.0000 \
  fxx:@:5:1.0000 \
  f0a:@:5:0.33333 \
  f1a:@:5:0.33333 \
  f0ax:-darklines@-linewd@2:5:0.33333

all: 
	for t in ${TESTS}; do \
          make SPECS=$$t single; \
        done

clean::

SPECS := SPECS.IS.UNDEFINED
ifneq "/${SPECS}" "/SPECS.IS.UNDEFINED"

SPECFIELDS := ${subst :, ,${SPECS}}
IMG := ${word 1,${SPECFIELDS}}
OPT := ${subst @, ,${word 2,${SPECFIELDS}}}
MAG := ${word 3,${SPECFIELDS}}
RED := ${word 4,${SPECFIELDS}}

single: ${IMG}-raw-g.ppm ${IMG}-fin.crs ${IMG}-fin-g.ppm

%-fin.crs: %-n.pgm %-raw.crs ${PROGDIR}/${PROG} Makefile
	cat $*-raw.crs \
          | egrep -v -e '^ *([#]|$$)' \
          | ${PROGDIR}/${PROG} \
              ${OPT} \
              -radius 7 \
              -maxdisp 10 -maxrot 8 \
              $*-n.pgm \
          > $*-fin.crs

%-raw-g.ppm: %-n.pgm %-raw.drw Makefile
	cat $*-n.pgm \
          | pnmscale ${MAG} \
          | convert PGM:- -fill yellow -draw '@'$*-raw.drw PPM:- \
          | pnmscale ${RED} \
          > $*-raw-g.ppm
	${SHOW} $*-raw-g.ppm $*-n.pgm

%-fin-g.ppm: %-n.pgm %-fin.drw Makefile
	cat $*-n.pgm \
          | pnmscale ${MAG} \
          | convert PGM:- -fill yellow -draw '@'$*-fin.drw PPM:- \
          | pnmscale ${RED} \
          > $*-fin-g.ppm
	${SHOW} $*-fin-g.ppm $*-n.pgm

%.drw: %.crs draw-crosses Makefile
	draw-crosses -v scale=${MAG} -v radius=${MAG} < $*.crs > $*.drw
        
endif
