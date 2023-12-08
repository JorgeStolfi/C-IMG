# Last edited on 2013-02-12 04:47:10 by stolfilocal

PROG := pnmwiggle
PROGDIR := ..

TESTS := \
  blk-1.pgm \
  ula-1.ppm
  
all:
	for file in ${TESTS}; do \
          do-one-test.sh $$file ; \
	done
  

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

clean::
	rm -f out/*.pgm out/*.ppm
        

