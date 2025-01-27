# Last edited on 2025-01-22 20:20:30 by stolfi

PROG := fni_view

PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# ${HTNAME}:${TXNAME}:${SCALE}

BUG_TESTS := \
  7-G-0256x0192.fni:200.000:NONE:NONE

ALL_TESTS := \
  ${BUG_TESTS} \
   \
  3-Z-0065x0049.fni:030.000:NONE:NONE \
  3-Z-0065x0049.fni:030.000:01-0065x0049.ppm:NONE \
  3-Z-0065x0049.fni:030.000:02-0065x0049.pgm:NONE \
  3-Z-0065x0049.fni:030.000:03-0064x0048.ppm:NONE \
  3-Z-0065x0049.fni:030.000:04-0064x0048.pgm:NONE \
   \
  3-G-0064x0048.fni:200.000:03-0064x0048.ppm:NONE \
  3-G-0064x0048.fni:200.000:03-0064x0048.ppm:1 \
   \
  4-N-0064x0048.fni:030.000:03-0064x0048.ppm:NONE \
  4-N-0064x0048.fni:030.000:NONE:0/1/2 \
   \
  5-G-0064x0048.fni:200.000:NONE:NONE \
   \
  6-G-0064x0048.fni:200.000:NONE:NONE \
  

NON_TESTS :=

TESTS := ${BUG_TESTS}
  
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
TXCHANS   := ${word 4, ${TEST_FIELDS}}

HTFILE = data/test-${HTNAME}

TX_OPS := 
ifeq "/${TXNAME}" "/NONE"
  TXFILE := 
  TX_OPS := 
else
  TXFILE := data/texture-${TXNAME}
  TX_OPS := -txFile ${TXFILE}
endif

ifeq "/${TXCHANS}" "/NONE"
  TXCHANS_OPS := 
else
  TXCHANS_OPS := -txChannels ${subst /, ,${TXCHANS}}
endif

run-single: ${PROGDIR}/${PROG} ${HTFILE} ${TXFILE}
	${PROGDIR}/${PROG} \
            ${TX_OPS} ${TXCHANS_OPS} \
            -channel 0 \
            -scale ${SCALE} \
            ${HTFILE}

endif
# End ${TEST_ARGS} section
######################################################################

