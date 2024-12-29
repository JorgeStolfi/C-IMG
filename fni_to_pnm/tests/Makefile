# Last edited on 2024-12-25 05:50:16 by stolfi
 
SHELL =	/bin/bash
       
PROG := fni_to_pnm

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# INNAME:OUTTAG:OUTEXT:OPTIONS

SOME_TESTS := \
  test-1-4:e:ppm:-channels¦0¦3¦7¦-center¦0¦-uniform

ALL_TESTS := \
  ramp-1up:a:pgm:-min¦0¦-max¦63¦-yAxis¦up \
  ramp-1dn:b:pgm:-min¦0¦-max¦63¦-yAxis¦down \
  ramp-1dn:c:pgm:-min¦0¦-max¦63 \
   \
  test-1-1:a:pgm:-channel¦0¦-min¦-16¦-max¦+16 \
  test-1-1:b:pgm:-channel¦0¦-center¦0¦-max¦+16 \
   \
  test-1-4:b:ppm:-channels¦0¦3¦7¦-min¦0¦0¦0¦-max¦1¦1¦1 \
  test-1-4:c:ppm:-channels¦0¦3¦7¦-min¦0¦-max¦1 \
  test-1-4:d:ppm:-channels¦0¦3¦7¦-min¦-2¦-1¦0¦-max¦1¦2¦3 \
  test-1-4:e:ppm:-channels¦0¦3¦7¦-center¦0¦-uniform

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
# Recursive "make" to create the test ${ZFUNC} with noise ${SIGMA}

ARGS := ${subst :, ,${TEST}}
NAME := ${word 1,${ARGS}}
TAG  := ${word 2,${ARGS}}
EXT  := ${word 3,${ARGS}}
OPTS := ${subst ¦, ,${word 4,${ARGS}}}

runtest: ${PROGDIR}/${PROG} in/${NAME}.fni Makefile
	@echo "NAME = ${NAME} OPTS = (${OPTS})"
	${PROGDIR}/${PROG} ${OPTS} \
          < in/${NAME}.fni > out/${NAME}-${TAG}.${EXT}
	display -title '%f' -filter Box -resize '<x600' out/${NAME}-${TAG}.${EXT}

# 
# End make ${ZFUNC}, ${SIGMA}
######################################################################
endif
