# Last edited on 2022-01-01 00:07:27 by stolfi

SHELL := /bin/bash

PROG := image_stitch
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# Test run arguments: {TESTNUM}:{VMAG}
# 
# where {TESTNUM} is test index and {VMAG} is the length scaling factor
# for error vector plotting.

TESTS := \
  4:03000.0 \
  5:00020.0

NONTESTS :=  \
  0:00020.0 \
  1:00020.0 \
  2:00020.0 \
  3:00020.0 \
  4:03000.0 \
  5:00020.0

.PHONY:: all-tests single-test

all: all-tests

all-tests:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST="$$tt" single-test ; \
        done
	
########################################################################
# Recursive make -- caller must define ${TEST}

TEST = TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

TESTARGS := ${subst :, ,${TEST}}
TESTNUM :=  ${word 1, ${TESTARGS}}
VMAG :=     ${word 2, ${TESTARGS}}

IPREFIX := in/test-${TESTNUM}
IPAIR_FILE := ${IPREFIX}-pairs.txt
PARM_FILE := ${IPREFIX}-parms.txt

OPREFIX := out/test-${TESTNUM}
OMAT1_FILE := ${OPREFIX}-1-matrix.txt
OMAT2_FILE := ${OPREFIX}-2-matrix.txt
OLIN1_FILE := ${OPREFIX}-1-outline.txt
OLIN2_FILE := ${OPREFIX}-2-outline.txt
OPAIR_FILE := ${OPREFIX}-pairs.txt
OPLOT_FILE := ${OPREFIX}-f2-plot.txt

single-test: ${PROGDIR}/${PROG} ${IPAIR_FILE} out
	-rm -f ${OPREFIX}-*.txt
	cat ${IPAIR_FILE} \
          | egrep -v '^[#]' \
          | ${PROGDIR}/${PROG} \
              `cat ${PARM_FILE} | sed -e 's:[\#].*$$::g'` \
              -writeStitched 1 \
              -outPrefix ${OPREFIX}
	if [[ -s ${OPLOT_FILE} ]]; then ( show-goalf-plot.sh ${TESTNUM} ${OPLOT_FILE} & ) ; xmessage "Ok" ; fi
	if [[ -s ${OPAIR_FILE} ]]; then show-mapping.sh ${TESTNUM} ${OLIN1_FILE} ${OLIN2_FILE} ${OPAIR_FILE} ${VMAG} ; fi


endif
# end ${TEST} make
########################################################################

clean::
	-/bin/rm -f out/test-*.{txt,png}
