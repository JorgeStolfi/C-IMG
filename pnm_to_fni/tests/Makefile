# Last edited on 2024-12-25 05:40:12 by stolfi

SHELL =	/bin/bash

PROG := pnm_to_fni
        
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# INNAME:OUTTAG:OUTEXT:OPTIONS

ALL_TESTS := \
  ramp-3:a:pgm:-min¦0¦-max¦63¦-yAxis¦down \
  ramp-3:b:pgm:-min¦0¦-max¦63¦-yAxis¦up\
  ramp-3:c:pgm:-min¦0¦-max¦63 \
   \
  erre-1:a:ppm:-min¦0¦-max¦-1¦-yAxis¦down \
  erre-1:b:ppm:-min¦0¦-max¦-1¦-yAxis¦up \
  erre-1:c:ppm:-min¦0¦-max¦-1 \
   \
  ramp-1:a:pgm:-min¦-16¦-max¦+16 \
  ramp-1:b:pgm:-center¦0¦-max¦+16 \
   \
  rbow-1:b:ppm:-min¦0¦0¦0¦-max¦1¦1¦1 \
  rbow-1:c:ppm:-min¦0¦-max¦1 \
  rbow-1:d:ppm:-min¦-2¦-1¦0¦-max¦1¦2¦3 \
   \
  

TESTS := ${ALL_TESTS}

.PHONY:: runtest-all runtest

all: runtest-all

# ====================================================================
# Cleanup

clean::
	-/bin/rm -f out/*.fni

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

runtest: ${PROGDIR}/${PROG} in/${NAME}.${EXT} Makefile
	@echo "NAME = ${NAME} OPTS = (${OPTS})"
	${PROGDIR}/${PROG} ${OPTS} \
          < in/${NAME}.${EXT} > out/${NAME}-${TAG}.fni
	( display -title '%f' -filter Box -resize '<800x600' in/${NAME}.${EXT} & ) ; \
          fni_view out/${NAME}-${TAG}.fni

# 
# End make ${ZFUNC}, ${SIGMA}
######################################################################
endif
