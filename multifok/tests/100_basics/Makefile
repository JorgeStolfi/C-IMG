# Last edited on 2017-10-26 01:35:49 by stolfilocal

PROG := multifok
PROGDIR := ../..

TESTS := \
  test3

.PHONY:: run-all-tests run-single debug-single show-single clean-single

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: run-all-tests

run-all-tests:
	for test in ${TESTS} ; do \
          ${MAKE} TEST_ARGS="$${test}" debug-single run-single ; \
        done
        
######################################################################
# For recursive makes - caller must define ${TEST_ARGS}
ifneq "/${TEST_ARGS}" "/"

TEST_ARGS_WDS := ${subst :, ,${TEST_ARGS}}

TEST_NAME  :=   ${word  1, ${TEST_ARGS_WDS}}

debug-single: 
	@echo "TEST_NAME = ${TEST_NAME}"
        
IMAGE_LIST := in/${TEST_NAME}-images.txt

OUT_PREFIX := out/${TEST_NAME}

OUT_FINAL_PNG := ${OUT_PREFIX}-final.png

.PRECIOUS: ${OUT_FINAL_PNG}

run-single: clean-single ${OUT_FINAL_PNG} show-single

${OUT_FINAL_PNG}: ${PROGDIR}/${PROG} 
	${PROGDIR}/${PROG} \
            `cat ${IMAGE_LIST} | egrep -v -e '^[ ]*([\#]|$$)'` \
            -verbose -noise 0.02 -spread 3.0 \
            -outPrefix ${OUT_PREFIX}

show-single:
ifeq "/${SHOW}" "/YES"
	display ${PLOT_FILE} 
endif

clean-single: 
	rm -fv ${OUT_PREFIX}*.{png,txt,fni}
        
endif
# end ${TEST_ARGS} section
######################################################################

clean::

