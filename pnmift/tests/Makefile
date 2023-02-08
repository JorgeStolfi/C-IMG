# Last edited on 2023-02-08 01:13:37 by stolfi

PROG := pnmift

PROGDIR  := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  footb-1:ppm \
  ladys-2:pgm

NON_TESTS := \
  figs1-1:pgm \
  ladys-1:pgm \
  ladys-2:pgm \
  mrng3-1:pgm \
  unifm-1:pgm \
  tree3-1:pgm \
  tree3-2:pgm\
  footb-1:ppm \
  
  
.PHONY:: clean-all run-all clean-single run-single

all: run-all

clean:: clean-all

run-all: ; \
  for tst in ${TESTS}; do \
    ${MAKE} TEST=$${tst} run-single; \
  done
  
clean-all: ; \
  for tst in ${TESTS}; do \
    ${MAKE} TEST=$${tst} clean-single; \
  done

######################################################################
# Recursive section -- client must define ${TEST}
TEST := 
ifneq "/${TEST}" "/"

TEST_ARGS := ${subst :, ,${TEST}}
TNAME := ${word 1, ${TEST_ARGS}}
EXT := ${word 2, ${TEST_ARGS}}

IMGNAME = ${shell echo ${TNAME} | sed -e 's/[-][0-9]$$//'}

SCRIPT := data/${TNAME}.doit
INFILE := data/${IMGNAME}.${EXT}
OT_PREFIX := out/${TNAME}

SHOWIMG := display -filter box -resize '200%'
PSVIEW := evince
# PSVIEW := gv

run-single: ${OT_PREFIX}-costs.pgm

${OT_PREFIX}-costs.pgm: ${PROGDIR}/${PROG} ${SCRIPT} ${INFILE}
	@echo SCRIPT = ${SCRIPT} INFILE = ${INFILE}
	${SCRIPT}
	-if [[ -s ${OT_PREFIX}-costs.pgm ]]; then \
	  ${SHOWIMG} ${INFILE} ${OT_PREFIX}*.p?m ; \
        fi
	-if [[ -s ${OT_PREFIX}-trees.eps ]]; then \
          ${PSVIEW} ${OT_PREFIX}-trees.eps ; \
        fi

clean-single:
	rm -f \
          ${OT_PREFIX}-labels.pgm \
          ${OT_PREFIX}-preds.pgm \
          ${OT_PREFIX}-roots.pgm \
          ${OT_PREFIX}-spread.${EXT} \
          ${OT_PREFIX}-extract.${EXT} \
          ${OT_PREFIX}-trees.eps \
          ${OT_PREFIX}-boxes.txt \

endif
# End ${TEST} section
######################################################################
