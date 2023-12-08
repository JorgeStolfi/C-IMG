# Last edited on 2011-06-05 22:05:11 by stolfi

PROG := langev

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  0050!0001!amg!127!0310!0080!amg-a!0200
  
TESTS_OK := \
  0002!0001!sma!127!0006!0003!sma-a!2000
  
TESTS_BUG :=

ALL_TESTS := \
  0002!0001!sma!127!0006!0003!sma-a!2000 \
  0010!0001!amg!127!0320!0175!amg-a!0200

.PHONY:: single single-clean single-run single-profile

ACTIONS := clean run profile
# ACTIONS := profile

all:
	for t in ${TESTS} ; do \
	  for ac in ${ACTIONS} ; do \
            ${MAKE} TEST="$$t" ACTION=$$ac single ; \
          done ; \
        done

clean:

######################################################################
# Recursive make - must define ${TEST} ${ACTION}

ifneq "/${TEST}" "/"
ifneq "/${ACTION}" "/"

FIELDS := ${subst !, ,${TEST}}

FRAMES  := ${word 1,${FIELDS}}
STEP    := ${word 2,${FIELDS}}
INNAME  := ${word 3,${FIELDS}}
LEVEL   := ${word 4,${FIELDS}}
STARTX  := ${word 5,${FIELDS}}
STARTY  := ${word 6,${FIELDS}}
OTNAME  := ${word 7,${FIELDS}}
RESIZE  := ${word 8,${FIELDS}}

INPREFIX := in-${INNAME}
OTPREFIX := out/${OTNAME}

PGMIN := ${INPREFIX}.pgm

PPMOT := ${OTPREFIX}-000000.ppm

.PRECIOUS: ${PPMOT}

single: single-${ACTION}

single-clean:
	rm -f ${OTPREFIX}-[0-9][0-9][0-9][0-9][0-9][0-9].ppm

single-run: ${PPMOT}
	display \
          -title '%f' \
          -delay 2 \
          -filter Box \
          -resize '${RESIZE}%' \
            ${PGMIN} \
            ${OTPREFIX}-[0-9][0-9][0-9][0-9][0-9][0-9].ppm 

${PPMOT}: ${PROGDIR}/${PROG} ${PGMIN} Makefile
	@echo "${PPMIN} ${COLOR} ${NL} ${NR} ${NT} ${NB} -> ${PPMOT}"
	${PROGDIR}/${PROG} \
            -generations ${FRAMES} -step ${STEP} \
            -relief ${PGMIN} ${LEVEL} \
            -start ${STARTX} ${STARTY} \
            -prefix ${OTPREFIX}
        
single-profile: gmon.out
	gprof ${PROGDIR}/${PROG} \
          --flat-profile=${PROG}-gprof-flat.txt \
          --graph=${PROG}-gprof-graph.txt \
          -I .. \
          -I ../../../JSLIBS/libjs \
          -I ../../../JSLIBS/libjspbm \
          -I ../../../JSLIBS/libimg \
          -I ../../../JSLIBS/libgeo
        
endif
endif
# end ${TEST} ${ACTION}
######################################################################
