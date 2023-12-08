# Last edited on 2018-09-10 18:51:44 by stolfilocal

PROG := multifok_analyze
PROGDIR := ..

TESTS := \
  syngrd_01

NON_TESTS := \
  flower_01 \
  flower_01_mono \
  syntha_01

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

FRAME_PREFIX := data/${TEST_NAME}/frame_
FRAME_SUFFIX := .png
FRAME_PATTERN := ${FRAME_PREFIX}%05d${FRAME_SUFFIX}
FRAME_FILES := ${sort ${wildcard ${FRAME_PREFIX}*[0-9]${FRAME_SUFFIX}}}
FRAME_ID_LIST := ${subst ${FRAME_PREFIX},,${subst ${FRAME_SUFFIX},,${FRAME_FILES}}}

REF_IMAGE := data/${TEST_NAME}/reference.png

OUT_DIR := out/${TEST_NAME}

OUT_COEFFS := ${OUT_DIR}/coeff.txt

.PRECIOUS: ${OUT_COEFFS}

run-single: clean-single ${OUT_COEFFS} show-single

${OUT_COEFFS}: ${PROGDIR}/${PROG} plot_observations.sh 
	mkdir -p ${OUT_DIR}
	${PROGDIR}/${PROG} \
            -framePattern "${FRAME_PATTERN}" \
            ${addprefix -frameID ,${FRAME_ID_LIST}} \
            -imageFormat PNG \
            -reference ${REF_IMAGE} \
            -windowSize 3 \
            -verbose \
            -noise 0.01 \
            -threshold 0.05 \
            -outDir ${OUT_DIR}
	plot_observations.sh SHOW ${OUT_DIR}

show-single: ${OUT_COEFFS}
	for id in ${FRAME_ID_LIST} ; do \
          ifile="${FRAME_PREFIX}$${id}${FRAME_SUFFIX}" ; \
          ofile="${OUT_DIR}/simil_$${id}.png" ; \
          qfiles1=( "${OUT_DIR}/simil_$${id}_q"{00..02}".png" ) ; \
          qfiles2=( "${OUT_DIR}/simil_$${id}_q"{03..05}".png" ); \
          convert +append $${ifile} $${qfiles1[@]} .temp1.png ; \
          convert +append $${ofile} $${qfiles2[@]} .temp2.png ; \
          convert -append .temp1.png .temp2.png .temp.png ; \
	  display -title "frame $${id}" -filter Box -resize '200%' .temp.png ; \
        done

clean-single: 
	rm -fv ${OUT_DIR}/*.{png,txt,fni}
        
endif
# end ${TEST_ARGS} section
######################################################################

clean::



