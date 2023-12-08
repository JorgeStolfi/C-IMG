# Last edited on 2010-03-28 19:33:31 by stolfi

PROG := pgmwmask
PROGDIR := ..

PLOT := YES

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# {WX}/{WY}/{SHAPE}/{KIND}/{SX}/{SY}/{PWR}/{VSELF}/{TRM_ORD}

TESTS_SIGMA_10 := \
  071/071/O/G/08.0/08.0/0.0/1.000/2 \
  041/041/O/U/00.0/00.0/0.0/1.000/2 \
  051/051/O/P/01.0/01.0/1.0/0.000/2 \
  099/099/O/P/01.0/01.0/2.0/0.000/2

TESTS := ${TESTS_SIGMA_10}
  
TESTS_SOON :=

TESTS_DONE := \
   \
  033/021/O/U/08.0/04.0/0.0/1.000/0 \
  033/021/O/U/08.0/04.0/0.0/1.000/1 \
  033/021/O/U/08.0/04.0/0.0/1.000/2 \
  033/021/R/U/08.0/04.0/0.0/1.000/0 \
  033/021/R/U/08.0/04.0/0.0/1.000/1 \
  033/021/R/U/08.0/04.0/0.0/1.000/2 \
   \
  033/021/O/G/08.0/04.0/0.0/1.000/0 \
  033/021/O/G/08.0/04.0/0.0/1.000/1 \
  033/021/O/G/08.0/04.0/0.0/1.000/2 \
  033/021/R/G/08.0/04.0/0.0/1.000/0 \
  033/021/R/G/08.0/04.0/0.0/1.000/1 \
  033/021/R/G/08.0/04.0/0.0/1.000/2 \
   \
  033/021/O/P/08.0/04.0/2.0/1.000/0 \
  033/021/O/P/08.0/04.0/2.0/1.000/1 \
  033/021/O/P/08.0/04.0/2.0/1.000/2 \
  033/021/R/P/08.0/04.0/2.0/1.000/0 \
  033/021/R/P/08.0/04.0/2.0/1.000/1 \
  033/021/R/P/08.0/04.0/2.0/1.000/2 \
   \
  033/021/O/P/03.0/02.0/1.0/1.000/0 \
  033/021/O/P/03.0/02.0/1.0/1.000/1 \
  033/021/O/P/03.0/02.0/1.0/1.000/2 \
  033/021/R/P/03.0/02.0/1.0/1.000/0 \
  033/021/R/P/03.0/02.0/1.0/1.000/1 \
  033/021/R/P/03.0/02.0/1.0/1.000/2 \
   \
  015/009/O/U/00.0/00.0/0.0/0.500/2 \
  015/009/R/U/00.0/00.0/0.0/0.500/2

TESTS_LATER :=

.PHONY:: test-single clean-single

ACTIONS := test-single
# ACTIONS := clean-single test-single

all:
	for t in ${TESTS}; do \
          ${MAKE} TESTARGS=$$t ${ACTIONS} ; \
        done

clean::
	rm -f out/*.pgm out/*.txt

######################################################################
# Submake for single test -- caller must define ${TESTARGS}

TESTARGS := TESTARGS.IS.UNDEFINED

ifneq "/${TESTARGS}" "/TESTARGS.IS.UNDEFINED"

FIELDS :=   ${subst /, ,${TESTARGS}}
WX :=         ${word 1,${FIELDS}}
WY :=         ${word 2,${FIELDS}}
WIN_TAG :=    ${word 3,${FIELDS}}
WHT_TAG :=    ${word 4,${FIELDS}}
SX :=         ${word 5,${FIELDS}}
SY :=         ${word 6,${FIELDS}}
PWR :=        ${word 7,${FIELDS}}
VSELF :=      ${word 8,${FIELDS}}
TRM_ORD :=    ${word 9,${FIELDS}}

SX_TAG :=    ${subst .,,${SX}}
SY_TAG :=    ${subst .,,${SY}}
PWR_TAG :=   ${subst .,,${PWR}}
VSELF_TAG := ${subst .,,${VSELF}}
TRM_TAG :=   t${TRM_ORD}

ifeq "/${WIN_TAG}" "/R"
  WIN_SPEC := rect ${WX} ${WY}
endif
ifeq "/${WIN_TAG}" "/O"
  WIN_SPEC := oval ${WX} ${WY}
endif

ifeq "/${WHT_TAG}" "/U"
  WHT_SPEC := uniform
endif
ifeq "/${WHT_TAG}" "/G"
  WHT_SPEC := gaussian ${SX} ${SY}
endif
ifeq "/${WHT_TAG}" "/P"
  WHT_SPEC := power ${SX} ${SY} ${PWR}
endif

ifeq "/${TRM_ORD}" "/0"
  TRM_SPEC := step
endif
ifeq "/${TRM_ORD}" "/1"
  TRM_SPEC := quadratic
endif
ifeq "/${TRM_ORD}" "/2"
  TRM_SPEC := biquadratic
endif

NAME := mask

OT_NAME := msk-${WHT_TAG}-d${SX_TAG}x${SY_TAG}-p${PWR_TAG}-s${VSELF_TAG}-${WIN_TAG}-w${WX}x${WY}-${TRM_TAG}

OT_IMAGE := out/${OT_NAME}.pgm
OT_TEXT := out/${OT_NAME}.txt

test-single: ${OT_IMAGE}

${OT_IMAGE}: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
              -window ${WIN_SPEC} \
              -weights ${WHT_SPEC} \
              -trimming ${TRM_SPEC} \
              -self ${VSELF} \
              -verbose \
              -maxval 65535 \
              -output pgm \
          > $@
	-display -filter box -resize '500%' $@
ifeq "/${PLOT}" "/YES"
	-plot_pgm_file.sh ${OT_IMAGE} 
endif

clean-single:
	rm -f ${OT_IMAGE}

endif
# End of ${TESTARGS}
######################################################################
