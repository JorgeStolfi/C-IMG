# Last edited on 2009-01-06 17:50:18 by stolfi

PROG := pnmadjust
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  ortho-200.adj

TESTS_NOT :=

all: ${TESTS}

clean::
	rm -f ${TESTS}

SHOW := display
HRULE := echo ' --------------------------------------------- '
        
%.adj: %.txt ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
          < $*.txt \
	  > $@
	cat $@
	@${HRULE}
