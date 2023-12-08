# Last edited on 2011-06-09 13:32:33 by stolfi

PROG := plan_extract
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  t1-50-50.txt \
  t2-50-50.txt

TESTS_NOT :=

all: ${TESTS}

clean::
	rm -f ${TESTS}

SHOW := display
HRULE := echo ' --------------------------------------------- '
        
out/t1-50-50.txt: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	  -input t1 640 480 \
          -minSize 160 120 \
          -midSize 224 168 \
          -maxSize 2560 1920
	cat $@
	@${HRULE}

out/t2-50-50.txt: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	  -input t2 640 480 \
          -take 635 467 2 1 \
          -minSize 640 480 \
          -maxSize 640 480 \
          -pad tbr \
          -maxPad 0.20 \
          -maxCut 0.20
	cat $@
	@${HRULE}

# -minSize 20 15
# -midSize 28 21
