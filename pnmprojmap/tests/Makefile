# Tests of the pnmprojmap tool
# Last edited on 2023-03-03 05:46:37 by stolfi

PROG := pnmprojmap
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DEBUG := -verbose -debug 140 195

TESTS := \
  ptpts:abi-1:ppm \
  ident:abi-1:ppm \
  rot90:flo-1:ppm \
  shr45:flo-1:ppm \
  rot45:flo-1:ppm \
  shf05:flo-1:ppm
  
TESTS_LATER_2 :=  \
  kmnpt:kmn-7:ppm \
  kmnpt:kmn-8:ppm \
  kmnpt:kmn-9:ppm  
 
TESTS_LATER := \
  ptpts:abi-1:ppm \
  ident:abi-1:ppm \
  persp:abi-1:ppm \
  rot90:flo-1:ppm \
  shr45:flo-1:ppm \
  rot45:flo-1:ppm \
  shf05:flo-1:ppm \

NONTESTS_GRAY := \
  ident:flo-1:pgm \
  shf05:flo-1:pgm \
  cubic:flo-1:pgm

NONTESTS_COLOR := \
  ptpts:abi-1:ppm \
  ident:abi-1:ppm \
  persp:abi-1:ppm \
  rot90:flo-1:ppm \
  shr45:flo-1:ppm \
  rot45:flo-1:ppm \
  shf05:flo-1:ppm

MAKEFILE := Makefile

.PHONY:: \
  run-tests single-test compare-images

all: clean run-tests

run-tests:
	for tt in ${TESTS}; do \
          ${MAKE} SPECS=$$tt single-test; \
        done

clean::
	rm -f out/*.p?m

SHOWIMG := display
# SHOWIMG := xv

########################################################################
# For recursive "make" of a single test. Caller must define ${SPECS}.

ifneq "/${SPECS}" "/"

SPECFIELDS := ${subst :, ,${SPECS}}
TST := ${word 1,${SPECFIELDS}}
IMG := ${word 2,${SPECFIELDS}}
EXT := ${word 3,${SPECFIELDS}}

INFILE := data/${IMG}.${EXT}
OTFILE := out/${IMG}-${TST}-out.${EXT}
RFFILE := out/${IMG}-${TST}-ref.${EXT}
DFFILE := out/${IMG}-${TST}-dif.${EXT}

single-test: ${OTFILE} ${RFFILE} ${DFFILE} ${INFILE}
	${SHOWIMG} ${OTFILE} ${RFFILE} ${DFFILE} ${INFILE}

${DFFILE}: ${OTFILE} ${RFFILE}
	pnmxarith -subtract -offset 0.5 ${OTFILE} ${RFFILE} > ${DFFILE}

#----------------------------------------------------------------------
# The identity map.

out/${IMG}-ident-ref.${EXT}: ${INFILE}
	cp -av ${INFILE} ${RFFILE}

out/${IMG}-ident-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -matrix \
              +1.000000000 00.000000000 00.000000000 \
              00.000000000 +1.000000000 00.000000000 \
              00.000000000 00.000000000 +1.000000000 \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

#----------------------------------------------------------------------
# Shift by half a pixel.

out/${IMG}-shf05-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -yAxis up \
            -matrix \
              +1.000000000 +0.500000000 +0.500000000 \
              00.000000000 +1.000000000 00.000000000 \
              00.000000000 00.000000000 +1.000000000 \
            -oCenter \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-shf05-ref.${EXT}: ${INFILE} ${MAKEFILE}
	cat ${INFILE} \
          > ${RFFILE}

#----------------------------------------------------------------------
# Shift by half a pixel, with bicubic interpolation.

out/${IMG}-cubic-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -yAxis up \
            -matrix \
              +1.000000000 +0.500000000 +0.500000000 \
              00.000000000 +1.000000000 00.000000000 \
              00.000000000 00.000000000 +1.000000000 \
            -oCenter \
            -interpolate 1 \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-cubic-ref.${EXT}: ${INFILE} ${MAKEFILE}
	cat ${INFILE} \
          > ${RFFILE}

#----------------------------------------------------------------------
# Rotation by 45 degrees.

out/${IMG}-rot45-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -yAxis up \
            -matrix \
              +1.000000000 00.000000000 00.000000000 \
              00.000000000 +0.707106781 +0.707106781 \
              00.000000000 -0.707106781 +0.707106781 \
            -oCenter \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-rot45-ref.${EXT}: ${INFILE} ${MAKEFILE}
	cat ${INFILE} \
          | pnmcut 1 1 399 399 \
          | pnmrotate 45 \
          | pnmcut 83 83 400 400 \
          > ${RFFILE}

#----------------------------------------------------------------------
# Rotation by 90 degrees.

out/${IMG}-rot90-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -matrix \
              +1.000000000 00.000000000 00.000000000 \
              00.000000000 00.000000000 -1.000000000 \
              00.000000000 +1.000000000 00.000000000 \
            -oCenter \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-rot90-ref.${EXT}: ${INFILE}
	cat ${INFILE} \
          | pamflip -r90 \
          > ${RFFILE}

#----------------------------------------------------------------------
# Shearing by 45 degrees.

out/${IMG}-shr45-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
            -yAxis up \
            -iCenter \
            -matrix \
              +1.000000000   00.000000000   00.000000000 \
              00.000000000   +1.000000000   00.000000000 \
              00.000000000   +1.000000000   +1.000000000 \
            -oCenter \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-shr45-ref.${EXT}: ${INFILE}
	cat ${INFILE} \
          | pnmshear -45 \
          | pnmcut 200 0 400 400 \
          > ${RFFILE}

#----------------------------------------------------------------------
# Perspective view given by matrix.

out/${IMG}-persp-ref.${EXT}: ${INFILE}
	cp -av ${INFILE} ${RFFILE}

out/${IMG}-persp-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -matrix \
              +1.500000000   00.000000000   00.000000000 \
              00.000000000   +0.900000000   00.000000000 \
              +0.002777777   00.000000000   +0.900000000 \
            -oCenter \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

#----------------------------------------------------------------------
# Perspective view given by 4 point pairs. 

out/${IMG}-ptpts-ref.${EXT}: ${INFILE}
	cp -av ${INFILE} ${RFFILE}

out/${IMG}-ptpts-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -points \
              010 010 \
              270 010 \
              010 382 \
              270 382 \
              \
              260 040 \
              260 280 \
              020 020 \
              020 300 \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

#----------------------------------------------------------------------

# Perspective view by 4 point pairs, from data/${IMG}.pts --> data/${IMG}-ref.pts. 

PTFILE := data/${IMG}.pts
PRFILE := data/${IMG}-ref.pts

out/${IMG}-kmnpt-ref.${EXT}: ${INFILE}
	cp -av ${INFILE} ${RFFILE}

out/${IMG}-kmnpt-out.${EXT}: ${INFILE} ${PTFILE} ${PRFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -yAxis down \
            -points \
              `cat ${PTFILE} | head -n 4 | sed -e 's:[\#].*$$::'` \
              `cat ${PRFILE} | head -n 4 | sed -e 's:[\#].*$$::'` \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

endif
# End  ${SPECS} section.
########################################################################

