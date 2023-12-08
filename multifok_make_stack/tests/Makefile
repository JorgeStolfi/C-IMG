# Last edited on 2018-07-25 15:16:38 by stolfilocal

PROG := multifok_make_stack
PROGDIR := ..

TESTS := \
  spikes_01 \
  syntha_01

MORE_TESTS := \
  syngrd_01

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

MAIN_IMAGE := data/${TEST_NAME}/main.png
HEIGHT_MAP := data/${TEST_NAME}/height.png

OUT_DIR := out/${TEST_NAME}
FRAME_PREFIX := ${OUT_DIR}/frame_
FRAME_PATTERN := ${FRAME_PREFIX}%05d

OUT_FRAME0 := ${FRAME_PREFIX}00000.png

.PRECIOUS: 

run-single: clean-single ${OUT_FRAME0} show-single

MIN_HEIGHT := -0.1000
HEIGHT_STEP := 0.1000
NUM_FRAMES := 13

${OUT_FRAME0}: ${PROGDIR}/${PROG}
	mkdir -pv ${OUT_DIR}
	${PROGDIR}/${PROG} \
            -imageFormat PNG \
            -image ${MAIN_IMAGE} \
            -heightMap ${HEIGHT_MAP} \
            -numFrames ${NUM_FRAMES} \
            -minHeight ${MIN_HEIGHT} \
            -heightStep ${HEIGHT_STEP} \
            -blurFactor 3.000 \
            -framePattern "${FRAME_PATTERN}" \
            -outDir ${OUT_DIR} \
            -verbose

show-single: ${OUT_FRAMES}
	for nid in {0..${NUM_FRAMES}} ; do \
	  id=`printf "%05d" "$${nid}"` ; \
          ffile="${FRAME_PREFIX}$${id}.png" ; \
          mfile="${FRAME_PREFIX}$${id}_mask.png" ; \
          echo "displaying $${ffile}..." ; \
          if [[ -s $${ffile} ]]; then \
            convert +append $${ffile} $${mfile} .temp.png ; \
	    display -title "frame $${id}" -filter Box -resize '400%' .temp.png ; \
          fi ; \
        done

clean-single: 
	rm -fv ${OUT_DIR}/*.{png,txt,fni}
        
endif
# end ${TEST_ARGS} section
######################################################################

clean::



