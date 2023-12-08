# Last edited on 2010-07-03 01:37:37 by stolfi

PROG := pnmadjust
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  air-canyon-1
  
NONTESTS := \
  2000-talossa-03 \
  eijistolfi-1 \
  jaburu-1

.PHONY:: single single-clean single-img

all:
	for t in ${TESTS} ; do \
          ${MAKE} TEST=$$t single ; \
        done

clean::
	for t in ${TESTS} ; do \
          ${MAKE} TEST=$$t single-clean ; \
        done

SHOW := display
HRULE := echo ' --------------------------------------------- '

######################################################################
# Recursive make - must define ${TEST}

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

PPMRAW := data/${TEST}.ppm.gz
JPGOLD := data/${TEST}-o.jpg

PPMCOR := out/${TEST}-c.ppm.gz
JPGCOR := out/${TEST}-c.jpg

single: single-clean single-img
	@${HRULE}

single-clean:
	rm -f ${PPMCOR} ${JPGCOR}

single-img: ${PPMCOR} ${JPGCOR} ${JPGOLD}
	${SHOW} ${PPMCOR} ${PPMRAW} ${JPGCOR} ${JPGOLD} 

ARGS := ${shell cat data/${TEST}.parms | sed -e 's/[\#].*//g' }

${PPMCOR}: ${PROGDIR}/${PROG}
	cat ${PPMRAW} \
          | gunzip \
          | ${PROGDIR}/${PROG} ${ARGS} \
          | gzip \
          > ${PPMCOR}

%.jpg: %.ppm.gz
	cat $< \
          | gunzip \
          | pnmgamma 2.2 \
          | convert -quality 98 PPM:- $*.jpg
        
endif
# end ${TEST}
######################################################################
