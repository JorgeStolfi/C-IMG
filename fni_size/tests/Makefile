# Last edited on 2023-03-29 13:21:28 by stolfi

PROG := fni_size

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  test-1:002/065/049 \
  test-2:000/000/000

.PHONY:: \
  runalltests runtest

# all:
all: runalltests

# ====================================================================
# Cleanup

clean::

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
NAME := ${word 1,${TESTFIELDS}}
REF  := ${word 2,${TESTFIELDS}}

IN_FILE := in/${NAME}.fni

runtest: ${PROGDIR}/${PROG} ${IN_FILE} Makefile
	@echo "${NAME}: ${REF}"
	@printf "  expected:    %s\n" "${REF}"
	@printf "  formatted:   %s\n" "`${PROGDIR}/${PROG} -format '%03d/%03d/%03d' < ${IN_FILE}`"
	@printf "  unformatted: %s\n" "`${PROGDIR}/${PROG} ${IN_FILE}`"

# 
# End make ${TESTARGS}
######################################################################
endif
