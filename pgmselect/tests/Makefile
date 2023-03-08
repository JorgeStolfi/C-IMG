# Last edited on 2023-03-03 14:55:14 by stolfi

PROG := pgmselect
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: clean do-test

TNAME := tt

INFILE := data/${TNAME}.pgm
OUTFILE := out/${TNAME}-sel.pgm

do-test: ${OUTFILE}

${OUTFILE}: ${PROGDIR}/${PROG} ${INFILE} Makefile
	${PROGDIR}/${PROG} \
          -range 100 150 \
          < ${INFILE} \
          > ${OUTFILE}
	if [[ -s ${OUTFILE} ]]; then \
          display ${INFILE} ${OUTFILE} ; \
        fi

clean::

