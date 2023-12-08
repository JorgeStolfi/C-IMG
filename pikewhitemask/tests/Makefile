# Last edited on 2023-02-04 07:07:24 by stolfi

PROG := pikewhitemask

PROGDIR  := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  t_000:501:1.0000:1.0000

NON_TESTS :=
  
.PHONY:: run-all run-single

all: run-all

clean:: 
	( cd out && rm -f *.ppm *.pgm )

run-all: ; \
  for tst in ${TESTS}; do \
    ${MAKE} TEST_ARGS=$${tst} run-single; \
  done

######################################################################
# Recursive section -- client must define ${TEST}
TEST_ARGS := 
ifneq "/${TEST_ARGS}" "/"

TEST_FIELDS := ${subst :, ,${TEST_ARGS}}
IMGNAME   := ${word 1, ${TEST_FIELDS}}
EQ_COL    := ${word 2, ${TEST_FIELDS}}
EQ_GAIN0  := ${word 3, ${TEST_FIELDS}}
EQ_GAIN1  := ${word 4, ${TEST_FIELDS}}

SCRIPT := data/${IMGNAME}.doit
INFILE := data/${IMGNAME}.raw
OTFILE := out/${IMGNAME}_mask.pgm

SHOWIMG := display
PSVIEW := evince
# PSVIEW := gv

run-single: ${OTFILE}

${OTFILE}: ${PROGDIR}/${PROG} ${INFILE}
	@echo INFILE = ${INFILE}
	${PROGDIR}/${PROG} \
            -equalize ${EQ_COL} ${EQ_GAIN0} ${EQ_GAIN1} \
            -verbose \
          < ${INFILE} \
          | pnmnoraw \
          > ${OTFILE}
	-${SHOWIMG} ${OTFILE}

endif
# End ${TEST_ARGS} section
######################################################################
