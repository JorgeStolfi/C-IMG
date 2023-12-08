# Last edited on 2011-11-21 09:04:32 by stolfi

PROG := stereo_adjust
PROGDIR := ..

all: do-test

do-test: ${PROGDIR}/${PROG} in/dados-1.txt
	cat in/dados-1.txt \
          | egrep -v '^[#]' \
          | ${PROG}  \
          > out/dados-1.txt

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

