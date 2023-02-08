# Last edited on 2012-02-25 09:55:03 by stolfilocal

PROG := ppmmarpap
PROGDIR := ..

PLOT := YES

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# {NX}:{NY}:{NW}:{ND}

TESTS := \
  1480:1050:064:016

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

FIELDS :=   ${subst :, ,${TESTARGS}}
NX :=         ${word 1,${FIELDS}}
NY :=         ${word 2,${FIELDS}}
NW :=         ${word 3,${FIELDS}}
ND :=         ${word 4,${FIELDS}}

NAME := mask

OT_NAME := pap

OT_IMAGE := out/${OT_NAME}.ppm

test-single: ${OT_IMAGE}

DATE := ${shell yyyy-mm-dd-hhmmss}

${OT_IMAGE}: ${PROGDIR}/${PROG} Makefile
	time ${PROGDIR}/${PROG} \
              -size ${NX} ${NY} \
              -whorls ${NW} \
              -detailSize ${ND} \
              -verbose \
              -maxval 65535 \
          > $@
	cp -av ${PROGDIR}/${PROG}.c out/${DATE}-${PROG}.c
	-convert $@ out/${DATE}-${OT_NAME}.png
	-display /tmp/init.p*m $@

clean-single:
	rm -f ${OT_IMAGE}

endif
# End of ${TESTARGS}
######################################################################
