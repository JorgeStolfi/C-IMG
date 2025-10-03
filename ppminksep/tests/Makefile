# Last edited on 2025-08-07 21:35:04 by stolfi

PROG := ppminksep

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  f1r-2-tian:a \
  f1r-3-big2:c
  
NONTESTS := \
  f81v-1-nym:a \
  f81v-2-brn:a \
  f81v-3-grn:b \
  f81v-4-txt:a \
  \
  f1r-1-oldy:a \
  f1r-2-tian:a \
  f1r-3-big2:c \
  \
  f4r-1-rot:d \
  f4r-2-vel:a \
  
TESTSBUG :=

TESTSOK :=

.PHONY:: single single-clean single-run

all:
	for t in ${TESTS} ; do \
          ${MAKE} TEST=$$t single ; \
        done

clean::
	for t in ${TESTS} ; do \
          ${MAKE} TEST=$$t single-clean ; \
        done

######################################################################
# Recursive make - must define ${TEST}

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

TARGS := ${subst :, , ${TEST}}
IMGNAME :=    ${word 1, ${TARGS}}
PARMTAG := ${word 2, ${TARGS}}

IN_PPM := data/${IMGNAME}.ppm
PARMS := data/${IMGNAME}-${PARMTAG}.parms

OUT_PREFIX := out/${IMGNAME}-${PARMTAG}
PPMBG := ${OUT_PREFIX}-bg.ppm
PPMCA := ${OUT_PREFIX}-ca.ppm
PPMCB := ${OUT_PREFIX}-cb.ppm
PPMSH := ${OUT_PREFIX}-sh.ppm
PPMEX := ${OUT_PREFIX}-ex.ppm
PPMXX := ${OUT_PREFIX}-xx.ppm

PPMOUTS := ${PPMBG} ${PPMCA} ${PPMCB} ${PPMSH} ${PPMEX} ${PPMXX} ${TIFXX} 

single: single-run
# single: single-clean single-run

single-clean:
	rm -f ${PPMOUTS}

single-run: single-show ${PPMBG}

single-show: ${PARMS} ${IN_PPM} Makefile
	view_ppm_colors.sh \
            ${IN_PPM} \
            ${PARMS}

${PPMBG}: ${PROGDIR}/${PROG} ${IN_PPM} ${PARMS} Makefile
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
              `cat ${PARMS} | sed -e 's/[#].*//g'` \
              -outPrefix ${OUT_PREFIX}
	if [[ -s ${PPMBG} ]]; then display -title '%f' ${PPMOUTS} ${IN_PPM} ; fi
        
endif
# end ${IMGNAME}
######################################################################
