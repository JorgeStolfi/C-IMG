# Last edited on 2023-03-03 15:22:21 by stolfi

PROG = make_bullseye
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
  
OUTPREFIX = \
  out/test-image

OUTFILE := ${OUTPREFIX}.pgm

all: ${OUTPREFIX}

${OUTPREFIX}: ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} -size 512 512 > ${OUTFILE}
	if [[ -s ${OUTFILE} ]]; then \
	  display -filter Box -resize '200%' ${OUTFILE} ; \
        fi

