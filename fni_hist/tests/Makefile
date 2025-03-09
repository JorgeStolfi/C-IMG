# Last edited on 2025-02-28 18:05:02 by stolfi

PROG := fni_hist.sh

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# NAME:CHANNEL:LOGSCALE:STEP:VMIN:VMAX

SOME_TESTS := \
  test-6:0:0:0.005:NONE:NONE

ALL_TESTS := \
  test-1:0:0:0.100:NONE:NONE \
  test-1:1:0:0.100:NONE:NONE \
  test-1:2:0:0.100:NONE:NONE \
  \
  test-3:0:0:0.250:NONE:NONE \
  test-3:1:0:0.250:NONE:NONE \
  test-3:1:0:0.250:-20.0:+3.7 \
  \
  test-4:0:0:0.250:NONE:NONE \
  test-4:0:1:0.250:NONE:NONE \
  \
  test-5:0:0:0.250:NONE:NONE \
  test-6:0:0:0.005:NONE:NONE
  
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

TESTARGS := TESTARGS.IS.UNDEFINED
ifneq "/${TESTARGS}" "/TESTARGS.IS.UNDEFINED"
#######################################################################
# Recursive "make" to run the test ${TESTARGS}

TESTFIELDS  := ${subst  :, ,${TESTARGS}}
NAME     := ${word 1,${TESTFIELDS}}
CHANNEL  := ${word 2,${TESTFIELDS}}
LOGSCALE := ${word 3,${TESTFIELDS}}
STEP     := ${word 4,${TESTFIELDS}}
VMIN     := ${word 5,${TESTFIELDS}}
VMAX     := ${word 6,${TESTFIELDS}}

ifneq "/${VMIN}" "/NONE"
  RANGE_OPT := -range ${VMIN} ${VMAX}
else
  RANGE_OPT := 
endif

ifeq "/${LOGSCALE}" "/1"
  LOGSCALE_OP := -logScale
else
  LOGSCALE_OP := 
endif

IN_FNI := in/${NAME}.fni
OUT_EPS := out/${NAME}-${CHANNEL}.eps

PSVIEW := evince

runtest: ${PROGDIR}/${PROG} ${IN_FNI}
	${PROGDIR}/${PROG} \
            -channel ${CHANNEL} \
            -step ${STEP} ${RANGE_OPT} \
            ${LOGSCALE_OP} \
            -title "${NAME} ch${CHANNEL} log${LOGSCALE}" \
          < ${IN_FNI} \
          > ${OUT_EPS}
	${PSVIEW} ${OUT_EPS}

# 
# End make ${TESTARGS}
######################################################################
endif
