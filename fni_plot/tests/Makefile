# Makefile for generic C program
# Last edited on 2025-02-16 22:46:34 by stolfi

PROG := fni_plot.sh

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  test-3:1 \
  test-6:0

ALL_TESTS := \
  test-3:0 \
  test-3:1 \
  test-4:0 \
  test-5:0 \
  test-6:0

.PHONY:: \
  runalltests runtest

# all:
all: runalltests

# ====================================================================
# Cleanup

clean::
	-/bin/rm -f out/*.png

# ====================================================================
# Running tests

runalltests: 
	mkdir -p out
	for t in ${TESTS} ; do \
          ${MAKE} TESTARGS=$$t runtest ; \
        done

TESTARGS := TESTARGS.IS.UNDEFINED
ifneq "/${TESTARGS}" "/TESTARGS.IS.UNDEFINED"
#######################################################################
# Recursive "make" to run the test ${TESTARGS}

TESTFIELDS  := ${subst  :, ,${TESTARGS}}
NAME    := ${word 1,${TESTFIELDS}}
CHANNEL := ${word 2,${TESTFIELDS}}

IN_FNI := in/${NAME}.fni
OUT_PNG := out/${NAME}-${CHANNEL}.eps

runtest: ${PROGDIR}/${PROG} ${IN_FNI}
	${PROGDIR}/${PROG} \
            -channel ${CHANNEL} \
            -size 600 600 \
            -title "${NAME}:${CHANNEL}" \
          < ${IN_FNI} \
          > ${OUT_PNG}
	display -title "${NAME}:${CHANNEL}" ${OUT_PNG}

# 
# End make ${TESTARGS}
######################################################################
endif
