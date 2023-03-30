# Last edited on 2023-03-29 13:27:30 by stolfi
 
SHELL =	/bin/bash
       
PROG := fni_to_png

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# INNAME:OUTTAG:OUTEXT:OPTIONS

SOME_TESTS := \
  test-1-4:e:-channels%0%3%7%-uniform%T \
  test-1-4:b:-channels%0%3%7%3%-min%-16%-8%-4%-2%-max%16%16%16%16

ALL_TESTS := \
  test-1-1:z: \
  test-1-1:y:-channel%0 \
  test-1-1:a:-channel%0%-min%0%-max%+16 \
  test-1-4:b:-channel%0%-min%-16%-max%+16 \
  test-1-4:c:-channel%1%-min%-16%-max%+16 \
  test-1-4:d:-channel%2%-min%-16%-max%+16 \
  test-1-4:e:-channel%3%-min%-16%-max%+16 \
  \
  test-1-4:a:-channels%0%3%-min%-16%-max%+16 \
  test-1-4:b:-channels%1%3%-min%-16%-max%+16 \
  \
  test-1-4:b:-channels%0%3%7%-min%0%0%0%-max%16%8%4 \
  test-1-4:c:-channels%0%3%7%-min%0%-max%16 \
  test-1-4:d:-channels%0%3%7%-min%-16%-8%0%-max%16%16%16 \
  test-1-4:e:-channels%0%3%7%-uniform%T \
  \
  test-1-4:b:-channels%0%3%7%3%-min%-16%-8%-4%-2%-max%16%16%16%16

#TESTS := ${ALL_TESTS}
TESTS := ${SOME_TESTS}

.PHONY:: runtest-all runtest

all: runtest-all

# ====================================================================
# Cleanup

clean::
	-/bin/rm -f out/*.png

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
# Recursive "make" to create the test ${ZFUNC} with noise ${SIGMA}

ARGS := ${subst :, ,${TEST}}
NAME := ${word 1,${ARGS}}
TAG  := ${word 2,${ARGS}}
OPTS := ${subst %, ,${word 3,${ARGS}}}

runtest: ${PROGDIR}/${PROG}
	@echo "NAME = ${NAME} OPTS = (${OPTS})"
	${PROGDIR}/${PROG} ${OPTS} \
          < in/${NAME}.fni > out/${NAME}-${TAG}.png
	display -title '%f' -filter Box -resize '<x600' out/${NAME}-${TAG}.png

# 
# End make ${ZFUNC}, ${SIGMA}
######################################################################
endif
