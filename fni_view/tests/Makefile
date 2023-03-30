# Last edited on 2023-03-29 13:27:14 by stolfi

PROG := fni_view

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# ${HTNAME}:${TXNAME}:${SCALE}

TESTS := \
  test-3-N.fni:030.000:NONE \
  test-0-F.ppm:001.000:test-3-N.fni \
  test-0-F.ppm:001.000:test-3-Z.fni \
  test-3-Z.fni:001.000:test-0-F.ppm \
  test-3-Z.fni:001.000:test-3-N.fni \

NON_TESTS :=
  
.PHONY:: run-all run-single

all: run-all

clean:: 

run-all: ; \
  for tst in ${TESTS}; do \
    ${MAKE} TEST_ARGS=$${tst} run-single; \
  done

######################################################################
# Recursive section -- client must define ${TEST}
TEST_ARGS := 
ifneq "/${TEST_ARGS}" "/"

TEST_FIELDS := ${subst :, ,${TEST_ARGS}}
HTNAME    := ${word 1, ${TEST_FIELDS}}
SCALE     := ${word 2, ${TEST_FIELDS}}
TXNAME    := ${word 3, ${TEST_FIELDS}}

HTFILE = data/${HTNAME}

TX_OPS := 
ifeq "/${TXNAME}" "/NONE"
  TXFILE := 
  TX_OPS := 
else
  TXFILE := data/${TXNAME}
  TX_OPS := -texture ${TXFILE}
endif

run-single: ${PROGDIR}/${PROG} ${HTFILE} ${TXFILE}
	${PROGDIR}/${PROG} \
            ${TX_OPS} \
            -channel 0 \
            -scale ${SCALE} \
            ${HTFILE}

endif
# End ${TEST_ARGS} section
######################################################################

