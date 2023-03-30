# Last edited on 2023-03-29 13:28:09 by stolfi

PROG := fni_hist.sh

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# NAME:CHANNEL:STEP:VMIN:VMAX

SOME_TESTS := \
  test-6:0:0.005:NONE:NONE

ALL_TESTS := \
  test-1:0:0.100:NONE:NONE \
  test-1:1:0.100:NONE:NONE \
  test-1:2:0.100:NONE:NONE \
  \
  test-3:0:0.250:NONE:NONE \
  test-3:1:0.250:NONE:NONE \
  test-3:1:0.250:-20.0:+3.7 \
  \
  test-4:0:0.250:NONE:NONE \
  test-5:0:0.250:NONE:NONE \
  test-6:0:0.005:NONE:NONE
  
TESTS := ${ALL_TESTS}

.PHONY:: \
  runalltests runtest

# all:
all: runalltests

# ====================================================================
# Cleanup

clean::
	-/bin/rm -f out/*.eps

# ====================================================================
# Running tests

runalltests: 
	for t in ${TESTS} ; do \
          ${MAKE} TESTARGS=$$t runtest ; \
        done
]
TESTARGS := TESTARGS.IS.UNDEFINED
ifneq "/${TESTARGS}" "/TESTARGS.IS.UNDEFINED"
#######################################################################
# Recursive "make" to run the test ${TESTARGS}

TESTFIELDS  := ${subst  :, ,${TESTARGS}}
NAME    := ${word 1,${TESTFIELDS}}
CHANNEL := ${word 2,${TESTFIELDS}}
STEP    := ${word 3,${TESTFIELDS}}
VMIN    := ${word 4,${TESTFIELDS}}
VMAX    := ${word 5,${TESTFIELDS}}

ifneq "/${VMIN}" "/NONE"
  RANGE_OPT := -range ${VMIN} ${VMAX}
else
  RANGE_OPT := 
endif

IN_FNI := in/${NAME}.fni
OUT_EPS := out/${NAME}-${CHANNEL}.eps

PSVIEW := evince

runtest: ${PROGDIR}/${PROG} ${IN_FNI}
	${PROGDIR}/${PROG} \
            -channel ${CHANNEL} \
            -step ${STEP} ${RANGE_OPT} \
            -title "${NAME}:${CHANNEL}" \
          < ${IN_FNI} \
          > ${OUT_EPS}
	${PSVIEW} ${OUT_EPS}

# 
# End make ${TESTARGS}
######################################################################
endif
