# Last edited on 2023-03-07 17:38:12 by stolfi

PROG := pgmcolorize
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# FNAME:RAWZERO:STYLE:CYCLES

# Note: style 0 (unsigned and signed) ignores the CYCLES argument.

TESTS := \
  grayramp:000/00255:0:1 \
  \
  grayramp:000/00255:0:0 \
  grayramp:000/00255:0:1 \
  grayramp:000/00255:0:2 \
   \
  earthsea:000/00255:0:1 \
  earthsea:000/00255:0:2
  
RAMP_TESTS := \
  grayramp:000/00255:0:1 \
  \
  grayramp:127/00255:0:1 \
  \
  grayramp:127/00255:1:0 \
  grayramp:127/00255:1:1 \
  grayramp:127/00255:1:2 \
  \
  grayramp:127/00255:2:0 \
  grayramp:127/00255:2:1 \
  grayramp:127/00255:2:2 \

OTHER_TESTS := \
  grayramp:127/00255:2:0 \
  grayramp:127/00255:2:1 \
  grayramp:127/00255:2:2 \
  \
  anagalia:127/00255:2:0 \
  anagalia:127/00255:2:1 \
  anagalia:127/00255:2:2 \
  \
  anagalia:127/00255:0:1 \
  \
  anagalia:127/00255:1:1 \
  anagalia:127/00255:2:1 \
  \
  venustop:000/00255:0:1 \
  venustop:000/00255:0:1 \
  \
  earthsea:000/00255:0:1 \
  \
  earthvis:000/00255:0:1 \
  \
  mooncrat:000/00255:0:1

NONTESTS := \
  grayramp:000/00255:1:0 \
  grayramp:000/00255:1:1 \
  grayramp:000/00255:1:2
  
TESTSBUG :=

TESTSOK :=

.PHONY:: single single-clean single-run

all:
	for t in ${TESTS} ; do \
          ${MAKE} TEST=$$t single-test ; \
        done

clean::

######################################################################
# Recursive make - must define ${TEST}

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

FIELDS := ${subst :, ,${TEST}}

FNAME   := ${word 1,${FIELDS}}
RAWZERO := ${word 2,${FIELDS}}
STYLE   := ${word 3,${FIELDS}}
CYCLES  := ${word 4,${FIELDS}}

ZERO := ${subst /, / ,${RAWZERO}}

IN_NAME := in/${FNAME}-in
IN_PGM_FILE := ${IN_NAME}.pgm

OUT_NAME := out/${FNAME}-${subst /,-,${RAWZERO}}-${STYLE}-${CYCLES}
OUT_PPM_FILE := ${OUT_NAME}.ppm

single-test: ${OUT_PPM_FILE}

${OUT_PPM_FILE}: ${IN_PGM_FILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
              -inGamma 1.0 \
              -zero ${ZERO} \
              -cycles ${CYCLES} \
              -style ${STYLE}  \
              -outGamma 2.5 \
          < ${IN_PGM_FILE} \
          > ${OUT_PPM_FILE}
	if [[ -s ${OUT_PPM_FILE} ]]; then \
	  display ${OUT_PPM_FILE} ${IN_PGM_FILE} ; \
        fi

single-clean:
	rm -f ${OUT_PPM_FILE}

endif
# end ${TEST}
######################################################################
