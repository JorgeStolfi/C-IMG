# Last edited on 2023-03-07 21:02:06 by stolfi

PROG := pnmxpad

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  abi@10@20@30@40@-white \
  abi@10@20@30@40@-black \
  abi@10@20@30@40@-color_rgbi:1.0/0.5/0.2 \
  abi@10@20@30@40@-self
  
NONTESTS :=
  
TESTSBUG :=

TESTSOK :=

.PHONY:: single single-clean single-run

all:
	for t in ${TESTS} ; do \
          ${MAKE} TEST="$$t" single ; \
        done

clean::

######################################################################
# Recursive make - must define ${TEST}

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

FIELDS := ${subst @, ,${TEST}}

NAME := ${word 1,${FIELDS}}
NL := ${word 2,${FIELDS}}
NR := ${word 3,${FIELDS}}
NT := ${word 4,${FIELDS}}
NB := ${word 5,${FIELDS}}
COLOR := ${subst _, ,${word 6,${FIELDS}}}

PPMIN := in/${NAME}.ppm
PPMOT := out/${NAME}-bd.ppm

PGMIN := in/${NAME}.pgm
PGMOT := out/${NAME}-bd.pgm

PBMIN := in/${NAME}.pbm
PBMOT := out/${NAME}-bd.pbm

# single: single-run
single: single-clean single-run

single-clean:
	rm -f ${PPMOT} ${PGMOT} ${PBMOT}

single-run: ${PPMOT} ${PGMOT} ${PBMOT}
	display -title '%f' ${PPMOT} ${PPMIN} ${PGMOT} ${PGMIN} ${PBMOT} ${PBMIN}

${PPMOT}: ${PROGDIR}/${PROG} ${PPMIN} Makefile
	@echo "${PPMIN} ${COLOR} ${NL} ${NR} ${NT} ${NB} -> ${PPMOT}"
	cat ${PPMIN} \
          | ${PROGDIR}/${PROG} \
              ${COLOR} -left ${NL} -right ${NR} -top ${NT} -bottom ${NB} \
          > ${PPMOT}

${PGMOT}: ${PROGDIR}/${PROG} ${PGMIN} Makefile
	@echo "${PGMIN} ${COLOR} ${NL} ${NR} ${NT} ${NB} -> ${PGMOT}"
	cat ${PGMIN} \
          | ${PROGDIR}/${PROG} \
              ${COLOR} -left ${NL} -right ${NR} -top ${NT} -bottom ${NB} \
          > ${PGMOT}

${PBMOT}: ${PROGDIR}/${PROG} ${PBMIN} Makefile
	@echo "${PBMIN} ${COLOR} ${NL} ${NR} ${NT} ${NB} -> ${PBMOT}"
	${PROGDIR}/${PROG} \
              ${COLOR} -left ${NL} -right ${NR} -top ${NT} -bottom ${NB} \
              ${PBMIN} \
          > ${PBMOT}
        
endif
# end ${TEST}
######################################################################
