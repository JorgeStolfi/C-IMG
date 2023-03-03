# Tests of the pnmgtran tool
# Last edited on 2023-03-03 05:53:37 by stolfi

PROG := pnmradist
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Pixel to debug (indices always from TOP left corner):
DEBUG := -verbose -debug 0 0

# Pixel sensor dimensions (millimeters):
HX := 0.04
HY := 0.04 

TESTS := \
  t1:svga:ppm:0:+0.000001 \
  t1:hvga:ppm:0:-0.000001 \
  t1:svga:ppm:0:+0.010000 \
  t1:hvga:ppm:0:-0.010000

NONTESTS_GRAY := \
  t1:hvga:pgm:0:+0.004000 \
  t1:hvga:pgm:0:-0.004000 \
  t1:svga:pgm:0:-0.002000 \
  t1:svga:pgm:0:+0.002000 \
  t1:svga:ppm:0:+0.010000 \
  t1:hvga:ppm:0:-0.010000

NONTESTS_COLOR :=

MAKEFILE := Makefile

.PHONY:: \
  run-tests single-test compare-images

all: clean run-tests

run-tests:
	for tt in ${TESTS}; do \
          ${MAKE} SPECS=$$tt single-test; \
        done

clean::
	rm -f out/*.ppm out/*.pgm 

SHOW_IMG := display -title '%f'
# SHOW_IMG := xv

########################################################################
# For recursive "make" of a single test. Caller must define 
#   ${SPECS} := ${TAG}:${IMG_NAME}:${IMG_EXT}:${EXTEND}:${KAPPA}

ifneq "/${SPECS}" "/"

SPECFIELDS  := ${subst :, ,${SPECS}}
TAG         := ${word 1,${SPECFIELDS}}
IMG_NAME    := ${word 2,${SPECFIELDS}}
IMG_EXT     := ${word 3,${SPECFIELDS}}
EXTEND      := ${word 4,${SPECFIELDS}}
KAPPA       := ${word 5,${SPECFIELDS}}

# Input file name:
IN_FILE := in/${IMG_NAME}.${IMG_EXT}

# Filename prefix for run-specific files:
RUN_PREFIX := ${IMG_NAME}-${TAG}-x${EXTEND}-k${KAPPA}

# Reference image ${REF_FILE} (if any):
# REF_FILE := in/${RUN_PREFIX}-ref.${IMG_EXT}
REF_FILE := 

# Output files: (0) copy of input, (1) mapped, (2) mapped and un-mapped
OUT_FILE_0 := out/${RUN_PREFIX}-ot0.${IMG_EXT}
OUT_FILE_1 := out/${RUN_PREFIX}-ot1.${IMG_EXT}
OUT_FILE_2 := out/${RUN_PREFIX}-ot2.${IMG_EXT}

# Error image ${RER_FILE} = ${OUT_FILE_1} - ${REF_FILE}:
# RER_FILE := out/${RUN_PREFIX}-rer.${IMG_EXT}
RER_FILE := 

# Inversion error image ${IER_FILE} = ${OUT_FILE_2} - ${OUT_FILE_0}:
IER_FILE := out/${RUN_PREFIX}-ier.${IMG_EXT}

# The "-extend" option:
EXTEND_OPTION := ${subst 1,-extend,${subst 0,,${EXTEND}}}

# Negated value of {KAPPA}:
NEG_KAPPA := ${shell echo "/${KAPPA}" | tr '+-' '-+' | tr -d '/' }

single-test: \
          ${OUT_FILE_0} \
          ${OUT_FILE_1} \
          ${OUT_FILE_2} \
          ${RER_FILE} \
          ${IER_FILE}
	${SHOW_IMG} \
          ${OUT_FILE_0} \
          ${OUT_FILE_1} \
          ${OUT_FILE_2} \
          ${REF_FILE} \
          ${RER_FILE} \
          ${IER_FILE} \
          ${IN_FILE}

# Create the .PGM file from the .PPM file: 
%.pgm: %.ppm
	convert $< $@

${OUT_FILE_0}: ${IN_FILE}
	cp -p ${IN_FILE} ${OUT_FILE_0}

${OUT_FILE_1}: ${OUT_FILE_0} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -kappa ${KAPPA} \
            -pixelSize ${HX} ${HY} \
            ${EXTEND_OPTION} \
            -oCenter \
            ${DEBUG} \
            ${OUT_FILE_0} \
          > ${OUT_FILE_1}

${OUT_FILE_2}: ${OUT_FILE_1} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            -iCenter \
            -kappa ${NEG_KAPPA} \
            -pixelSize ${HX} ${HY} \
            ${EXTEND_OPTION} \
            -oCenter \
            ${DEBUG} \
            ${OUT_FILE_1} \
          > ${OUT_FILE_2}

${IER_FILE}: ${OUT_FILE_0} ${OUT_FILE_2}
	pnmxarith -subtract -offset 128 ${OUT_FILE_2} ${OUT_FILE_0} > ${IER_FILE}

########################################################################
# Section that compares ${OUT_FILE_1} with reference image ${REF_FILE}.
# Client must define ${REF_FILE} and${RER_FILE}.

ifneq "/" "/"
ifneq "/" "/"

${RER_FILE}: ${OUT_FILE_1} ${REF_FILE}
	pnmxarith -subtract -offset 128 ${OUT_FILE_1} ${REF_FILE} > ${RER_FILE}

endif 
endif

# End ${REF_FILE}, ${RER_FILE} section
########################################################################

endif
# End ${SPECS} section.
########################################################################
