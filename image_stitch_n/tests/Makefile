# Last edited on 2012-01-02 22:27:04 by stolfi

SHELL := /bin/bash

PROG := image_stitch_n
PROGDIR := ..

all: do-test

TEST := 2
INFILES := ${wildcard in/test-${TEST}-*-*-pairs.txt}
SHFILE := in/test-${TEST}-run.sh
OTFILE := out/test-${TEST}-matrix.txt
PLOTFILE := out/plot.txt

do-test: ${PROGDIR}/${PROG} ${INFILES} out
	-rm -f ${PLOTFILE}
	export PROGDIR=${PROGDIR} PROG=${PROG}; ${SHFILE}
	if [[ -s ${PLOTFILE} ]]; then show-goalf-plot.sh; fi

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

