# Last edited on 2023-02-08 11:36:10 by stolfi

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

LOCAL_PARMS := ${OUT_PREFIX}.parms

PPMOUT := ${PPMBG} ${PPMCA} ${PPMCB} ${PPMSH} ${PPMEX} ${PPMXX} ${TIFXX} 

single: single-run
# single: single-clean single-run

single-clean:
	rm -f ${PPMOUT} ${PPMOUTGZ}

single-run: ${PPMBG} single-show

${LOCAL_PARMS}: ${PARMS} Makefile
	cat ${PARMS} \
           | sed \
               -e 's/[#].*//g' \
           > ${LOCAL_PARMS}

single-show: ${LOCAL_PARMS} ${PPMBG} Makefile
	display -title '%f' ${PPMOUT} ${IN_PPM}
	./show_color_cloud.sh \
            `cat ${LOCAL_PARMS} | egrep '[-](inGamma|logScale)'` \
            ${IN_PPM} \
            ${OUT_PREFIX}

${PPMBG}: ${PROGDIR}/${PROG} ${IN_PPM} ${LOCAL_PARMS} Makefile
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
              `cat ${LOCAL_PARMS}` \
              -outPrefix ${OUT_PREFIX}
        
endif
# end ${IMGNAME}
######################################################################
