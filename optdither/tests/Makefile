# Last edited on 2023-03-17 22:16:57 by stolfi

PROG := optdither
PROGDIR := ..

.PHONY:: run-test

all: run-test

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

IN_FILE := in/pmat-knuth-08x08.dat

OUT_PREFIX := out/test

run-test: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
          -inFile ${IN_FILE} \
          -filterFreq 3.0 \
          -outPrefix ${OUT_PREFIX}
