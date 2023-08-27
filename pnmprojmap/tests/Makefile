# Tests of the pnmprojmap tool
# Last edited on 2023-08-26 23:57:43 by stolfi

PROG := pnmprojmap
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DEBUG := -verbose -debug 140 195

TESTS_BASIC := \
  ptpts:abi-1:ppm:1.000 \
  ident:abi-1:ppm:1.000 \
  rot90:flo-1:ppm:1.000 \
  shr45:flo-1:ppm:1.000 \
  rot45:flo-1:ppm:1.000 \
  shf05:flo-1:ppm:1.000
  
TEST_BASIC_SCALE := \
  ptpts:abi-1:ppm:0.250 \
  ident:abi-1:ppm:0.250 \
  rot90:flo-1:ppm:0.250 \
  shr45:flo-1:ppm:0.250 \
  rot45:flo-1:ppm:0.250 \
  shf05:flo-1:ppm:0.250 
  
TESTS_KMNPT := \
  kmnpt:kmn-7:ppm:1.000 \
  kmnpt:kmn-8:ppm:1.000 \
  kmnpt:kmn-9:ppm:1.000 \
   \
  kmnpt:kmn-7:ppm:1.500 \

NONTESTS_GRAY := \
  ident:flo-1:pgm:1.000 \
  shf05:flo-1:pgm:1.000 \
  cubic:flo-1:pgm:1.000

NONTESTS_COLOR := \
  ptpts:abi-1:ppm:1.000 \
  ident:abi-1:ppm:1.000 \
  persp:abi-1:ppm:1.000 \
  rot90:flo-1:ppm:1.000 \
  shr45:flo-1:ppm:1.000 \
  rot45:flo-1:ppm:1.000 \
  shf05:flo-1:ppm:1.000

# TESTS := ${TESTS_BASIC} ${TESTS_BASIC_SCALE} ${TESTS_KMNPT}
TESTS := ${TESTS_BASIC} ${TESTS_BASIC_SCALE}
# TESTS := ${TESTS_KMNPT}

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
TST :=   ${word 1,${SPECFIELDS}}
IMG :=   ${word 2,${SPECFIELDS}}
EXT :=   ${word 3,${SPECFIELDS}}
SCALE := ${word 4,${SPECFIELDS}}

SCALETAG := ${subst .,,${SCALE}}
INFILE := data/${IMG}.${EXT}
OTFILE := out/${IMG}-${TST}-${SCALETAG}-out.${EXT}
RFFILE := out/${IMG}-${TST}-${SCALETAG}-ref.${EXT}
DFFILE := out/${IMG}-${TST}-${SCALETAG}-dif.${EXT}

single-test: ${OTFILE} ${RFFILE} ${DFFILE} ${INFILE}
	${SHOWIMG} ${OTFILE} ${RFFILE} ${DFFILE} ${INFILE}

${DFFILE}: ${OTFILE} ${RFFILE}
	pnmxarith -subtract -offset 0.5 ${OTFILE} ${RFFILE} > ${DFFILE}
        
MAP_UNMAP := \
  -map P1 10 10 \
  -unmap Q1 10 10 \
  -map P2 20 20 \
  -unmap Q2 20 20 

#----------------------------------------------------------------------
# The identity map.

out/${IMG}-ident-${SCALETAG}-ref.${EXT}: ${INFILE}
	cat ${INFILE} \
          | pnmscale ${SCALE} \
          > ${RFFILE}

out/${IMG}-ident-${SCALETAG}-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -matrix \
              +1.000000000 00.000000000 00.000000000 \
              00.000000000 +1.000000000 00.000000000 \
              00.000000000 00.000000000 +1.000000000 \
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

#----------------------------------------------------------------------
# Shift by half a pixel.

out/${IMG}-shf05-${SCALETAG}-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -yAxis up \
            -matrix \
              +1.000000000 +0.500000000 +0.500000000 \
              00.000000000 +1.000000000 00.000000000 \
              00.000000000 00.000000000 +1.000000000 \
            -oCenter \
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-shf05-${SCALETAG}-ref.${EXT}: ${INFILE} ${MAKEFILE}
	cat ${INFILE} \
          | pnmscale ${SCALE} \
          > ${RFFILE}

#----------------------------------------------------------------------
# Shift by half a pixel, with bicubic interpolation.

out/${IMG}-cubic-${SCALETAG}-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -yAxis up \
            -matrix \
              +1.000000000 +0.500000000 +0.500000000 \
              00.000000000 +1.000000000 00.000000000 \
              00.000000000 00.000000000 +1.000000000 \
            -oCenter \
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            -interpolate 1 \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-cubic-${SCALETAG}-ref.${EXT}: ${INFILE} ${MAKEFILE}
	cat ${INFILE} \
          | pnmscale ${SCALE} \
          > ${RFFILE}

#----------------------------------------------------------------------
# Rotation by 45 degrees.

out/${IMG}-rot45-${SCALETAG}-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -yAxis up \
            -matrix \
              +1.000000000 00.000000000 00.000000000 \
              00.000000000 +0.707106781 +0.707106781 \
              00.000000000 -0.707106781 +0.707106781 \
            -oCenter \
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-rot45-${SCALETAG}-ref.${EXT}: ${INFILE} ${MAKEFILE}
	cat ${INFILE} \
          | pnmcut 1 1 399 399 \
          | pnmrotate 45 \
          | pnmcut 83 83 400 400 \
          | pnmscale ${SCALE} \
          > ${RFFILE}

#----------------------------------------------------------------------
# Rotation by 90 degrees.

out/${IMG}-rot90-${SCALETAG}-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -matrix \
              +1.000000000 00.000000000 00.000000000 \
              00.000000000 00.000000000 -1.000000000 \
              00.000000000 +1.000000000 00.000000000 \
            -oCenter \
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-rot90-${SCALETAG}-ref.${EXT}: ${INFILE}
	cat ${INFILE} \
          | pamflip -r90 \
          | pnmscale ${SCALE} \
          > ${RFFILE}

#----------------------------------------------------------------------
# Shearing by 45 degrees.

out/${IMG}-shr45-${SCALETAG}-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
            -yAxis up \
            -iCenter \
            -matrix \
              +1.000000000   00.000000000   00.000000000 \
              00.000000000   +1.000000000   00.000000000 \
              00.000000000   +1.000000000   +1.000000000 \
            -oCenter \
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-shr45-${SCALETAG}-ref.${EXT}: ${INFILE}
	cat ${INFILE} \
          | pnmshear -45 \
          | pnmcut 200 0 400 400 \
          | pnmscale ${SCALE} \
          > ${RFFILE}

#----------------------------------------------------------------------
# Perspective view given by matrix.

out/${IMG}-persp-${SCALETAG}-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -matrix \
              +1.500000000   00.000000000   00.000000000 \
              00.000000000   +0.900000000   00.000000000 \
              +0.002777777   00.000000000   +0.900000000 \
            -oCenter \
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-persp-${SCALETAG}-ref.${EXT}: ${INFILE}
	cat ${INFILE} \
          | pnmscale ${SCALE} \
          > ${RFFILE}

#----------------------------------------------------------------------
# Perspective view given by 4 point pairs. 

out/${IMG}-ptpts-${SCALETAG}-out.${EXT}: ${INFILE} ${PROGDIR}/${PROG}
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
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

out/${IMG}-ptpts-${SCALETAG}-ref.${EXT}: ${INFILE}
	cat ${INFILE} \
          | pnmscale ${SCALE} \
          > ${RFFILE}

#----------------------------------------------------------------------

# Perspective view by 4 point pairs, from data/${IMG}.pts --> data/${IMG}-ref.pts. 
# Assumes that the ref image is "kmn-9.ppm" no matter what the input image.

PTFILE := data/${IMG}.pts
PRFILE := data/${IMG}-ref.pts

out/${IMG}-kmnpt-${SCALETAG}-out.${EXT}: ${INFILE} ${PTFILE} ${PRFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -yAxis down \
            -points \
              `cat ${PTFILE} | head -n 4 | sed -e 's:[\#].*$$::'` \
              `cat ${PRFILE} | head -n 4 | sed -e 's:[\#].*$$::'` \
            -scale ${SCALE} \
            ${MAP_UNMAP} \
            ${DEBUG} \
            ${INFILE} \
          > ${OTFILE}

IRFILE := data/kmn-9.ppm

out/${IMG}-kmnpt-${SCALETAG}-ref.${EXT}: ${IRFILE}
	cat ${IRFILE} \
          | pnmscale ${SCALE} \
          > ${RFFILE}

endif
# End  ${SPECS} section.
########################################################################

