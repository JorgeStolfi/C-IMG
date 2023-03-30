# Last edited on 2023-03-29 18:31:09 by stolfi
 
SHELL =	/bin/bash
       
PROG := fni_to_pov

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# INNAME:OUTTAG:OUTEXT:OPTIONS

ALL_TESTS := \
  test-1:a: \
  test-1:b:-scaleZ%3 \
  test-2:a:

TESTS := ${ALL_TESTS}

.PHONY:: runtest-all runtest

all: runtest-all

# ====================================================================
# Cleanup

clean::
	-/bin/rm -f out/*.pgm out/*.ppm out/*.pbm 

# ====================================================================
# Running tests

runtest-all: 
	mkdir -p out
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt runtest ; \
        done

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"
#######################################################################
# Recursive "make" -- caller must define ${TEST}

ARGS := ${subst :, ,${TEST}}
NAME := ${word 1,${ARGS}}
TAG :=  ${word 2,${ARGS}}
OPTS := ${subst %, ,${word 3,${ARGS}}}

runtest: ${PROGDIR}/${PROG} in/${NAME}.fni Makefile
	@echo "NAME = ${NAME} OPTS = (${OPTS})"
	${PROGDIR}/${PROG} ${OPTS} \
          < in/${NAME}.fni \
          > out/${NAME}-${TAG}.inc
	( cd out && ../run_fni_povray.sh ../main.pov ${NAME}-${TAG} )	

# 
# End make ${ARGS}
######################################################################
endif
