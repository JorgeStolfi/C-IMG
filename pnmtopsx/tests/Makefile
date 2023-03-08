# Last edited on 2023-03-07 18:31:45 by stolfi

PROG := pnmtopsx
PROGDIR  := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

PICNAME := lady3

# PSVIEW := gv
PSVIEW := evince

TESTS := \
  pbm:00001: \
  pbm:00001:-scale@0.25 \
  pgm:00019: \
  pgm:00255: \
  ppm:00019: \
  ppm:00255: \
  ppm:00255:-turn \
  ppm:00255:-noturn \
  ppm:00255:-scale@0.50 \
  ppm:00255:-rle
  
TESTS_OK := \
  pbm:00001: \
  pbm:00001:-scale@0.25 \
   \
  pgm:00019: \
  pgm:00019:-scale@0.25 \
   \
  pgm:00019: \
  pgm:00255: \
  ppm:00019: \
  ppm:00255: \
  ppm:00255:-turn \
  ppm:00255:-noturn \
  ppm:00255:-scale@0.50 \
  ppm:00255:-rle
  
TESTS_DONT := \
  pgm:65535:@ \
  ppm:65535:@

ACTIONS := clean build-eps show-eps

all:
	for t in ${TESTS} ; do \
          ${MAKE} TESTARGS=$$t ${addprefix single-,${ACTIONS}} ; \
        done

clean::
	${MAKE} ACTIONS=clean  all

######################################################################
# Recursive make - caller must define ${TESTARGS}.

TESTARGS := TESTARGS.IS.UNDEFINED

ifneq "/${TESTARGS}" "/TESTARGS.IS.UNDEFINED"

TESTFIELDS := ${subst :, ,${TESTARGS}}
EXT :=    ${word 1,${TESTFIELDS}}
DEPTH :=  ${word 2,${TESTFIELDS}}
OPTS :=   ${word 3,${TESTFIELDS}}

OPTSTAG := o${subst .,,${subst @,-,${OPTS}}}

PNMNAME := ${PICNAME}-${DEPTH}
EPSNAME := ${PICNAME}-${DEPTH}-${EXT}-${OPTSTAG}

PNGFILE := in/${PICNAME}.png
PNMFILE := /tmp/${PNMNAME}.${EXT}
EPSFILE := out/${EPSNAME}.eps

single-build-eps: ${EPSFILE}

${EPSFILE}: ${PNMFILE} ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
            ${subst @, ,${OPTS}} \
            ${PNMFILE} \
          > ${EPSFILE}
        
single-show-eps: ${EPSFILE}
	${PSVIEW} ${EPSFILE}

/tmp/${PNMNAME}.pbm: ${PNGFILE}
	convert ${PNGFILE} PGM:- \
          | pgmnorm -bpercent 5 -wpercent 10 \
          | pgmtopbm -floyd \
          > $@
	pnmfile $@
	display $@ ${PNGFILE}

/tmp/${PNMNAME}.pgm: ${PNGFILE}
	convert ${PNGFILE} PGM:- \
          | pnmdepth ${DEPTH} \
          > $@
	pnmfile $@
	display $@

/tmp/${PNMNAME}.ppm: ${PNGFILE}
	convert ${PNGFILE} PPM:- \
          | pnmdepth ${DEPTH} \
          > $@
	pnmfile $@
	display $@

single-clean: 
	rm -f ${PNMFILE} ${EPSFILE}

endif
# end ${TESTARGS} section
######################################################################
       
