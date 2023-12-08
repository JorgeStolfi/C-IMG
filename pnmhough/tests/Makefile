# Last edited on 2011-06-05 23:55:29 by stolfi

PROG := pnmhough
PROGDIR := ..

TESTE := 6
EXT := pgm

.PHONY:: roda-teste

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: roda-teste

clean::
	rm out/*.pgm

roda-teste:
	${PROGDIR}/${PROG} teste-${TESTE}-in.${EXT} out/teste-${TESTE} 9
	display out/teste-${TESTE}*.p*m
