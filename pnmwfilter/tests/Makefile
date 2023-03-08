# Last edited on 2017-07-28 04:15:39 by jstolfi

PROG := pnmwfilter
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

SHOW := YES
GPROF := NO

# Format: "{IMG_NAME}:{MASK_NAME}:{FILTER}@{FILTER_ARGS}@.."

TEST_IMAGES_SOME_GRAY := \
  i0010:pgm \
  i1002:pgm \
  i2001:pgm

TEST_IMAGES_SOME_COLOR := \
  i3001:ppm \
  i3002:ppm \
  i3003:ppm
 
TEST_IMAGES_ALL_GRAY := \
  i0010:pgm \
  i0037:pgm \
  i0045:pgm \
  i1000:pgm \
  i1002:pgm \
  i1006:pgm \
  i1009:pgm \
  i2001:pgm
 
TEST_IMAGES_ALL_COLOR := \
  i3001:ppm \
  i3002:ppm \
  i3003:ppm \
  i3004:ppm \
  i3005:ppm \
  i3006:ppm \
  i3007:ppm \
  i3008:ppm \
  i3009:ppm
 
# TEST_IMAGES := ${TEST_IMAGES_ALL_GRAY}
TEST_IMAGES := ${TEST_IMAGES_SOME_GRAY}
# TEST_IMAGES := ${TEST_IMAGES_ALL_COLOR}

BASIC_FILTER_TESTS := \
  gauss-030x030:average \
  gauss-030x030:deviation \
  gauss-030x030:normalize@3.0 \
  gauss-030x030:rank \
  gauss-030x030:percentile@0.00 \
  gauss-030x030:percentile@0.00001 \
  gauss-030x030:percentile@0.25 \
  gauss-030x030:median \
  gauss-030x030:percentile@0.75 \
  gauss-030x030:percentile@0.99999 \
  gauss-030x030:percentile@1.00 \
  gauss-030x030:stretch@0.0001@0.9999 \
  gauss-030x030:stretch@0.25@0.75 \
  gauss-030x030:stretch@0.05@0.95

MASKED_FILTER_TESTS := \
  gauss-030x030:normalize@3.0@-mask@mk/m1000.pgm

SORT_FILTER_TESTS := \
  gauss-020x020:median

NONSORT_FILTER_TESTS := \
  gauss-030x030:average

NONSQUARE_KERNEL_TESTS := \
  gauss-040x020:normalize@3.0

KERNEL_SIZE_TESTS := \
  gauss-020x020:normalize@3.0 \
  gauss-030x030:normalize@3.0 \
  gauss-040x040:normalize@3.0

NOISE_LEVEL_TESTS := \
  gauss-030x030:stretch@0.05@0.95@-noise@000.0 \
  gauss-030x030:stretch@0.05@0.95@-noise@004.0 \
  gauss-030x030:stretch@0.05@0.95@-noise@064.0 \
  gauss-030x030:stretch@0.05@0.95@-noise@512.0

KERNEL_KIND_TESTS := \
  biqua-072x072:normalize@3.0 \
  gauss-072x072:normalize@3.0 \
  pwr10-072x072:normalize@3.0 \
  pwr20-072x072:normalize@3.0

# TESTS := ${NOISE_LEVEL_TESTS}
# TESTS := ${BASIC_FILTER_TESTS}
# TESTS := ${MASKED_FILTER_TESTS}
TESTS := ${KERNEL_KIND_TESTS}

ALL_TESTS := \
  ${BASIC_FILTER_TESTS} \
  ${SORT_FILTER_TESTS} \
  ${NONSORT_FILTER_TESTS} \
  ${NONSQUARE_KERNEL_TESTS} \
  ${KERNEL_SIZE_TESTS} \
  ${KERNEL_KIND_TESTS} \
  ${NOISE_LEVEL_TESTS}

.PHONY:: \
  test-set clean-set \
  test-single clean-single

ACTIONS := test-set
# ACTIONS := clean-set test-set

all:
	for t in ${TESTS}; do \
          ${MAKE} TESTARGS="$$t" ${ACTIONS} ; \
        done

######################################################################
# Submake for test on various images -- caller must define ${TESTARGS}

TESTARGS := 

ifneq "/${TESTARGS}" "/"

TESTFIELDS := ${subst :, ,${TESTARGS}}
MASK_NAME :=  ${word 1,${TESTFIELDS}}
RAW_SPECS :=  ${word 2,${TESTFIELDS}}

IN_DIR := in
WT_DIR := wt
OT_DIR := out

FILTER_SPECS := ${subst @, ,${RAW_SPECS}}
FILTER_NAME := ${word 1,${FILTER_SPECS}}

# Weight mask:
WT_IMAGE := ${WT_DIR}/${MASK_NAME}.pgm

# Tag for output file names:
TMP_XXX := ${subst mk/,,${subst .,,${subst @,-,${subst +,p,${subst -,m,${RAW_SPECS}}}}}}
FILTER_TAG := ${subst ${FILTER_NAME}-,,${TMP_XXX}}

OT_PREFIX := ${OT_DIR}/${FILTER_NAME}/${MASK_NAME}-${FILTER_TAG}
OT_STRIP_IMAGE := ${OT_PREFIX}-strip.png

test-set:
	rm -f ${OT_PREFIX}-i*.p?m
	for img in ${TEST_IMAGES} ; do \
          ${MAKE} TESTARGS="${TESTARGS}" IMGARGS="$$img" test-single ; \
        done
	convert ${OT_PREFIX}-i*.p?m +append ${OT_STRIP_IMAGE}
ifeq "/${SHOW}" "/YES"
	-if [ -s ${OT_IMAGE} ] ; then \
          display -title '%d/%f' ${OT_STRIP_IMAGE} ; \
        fi
endif  

clean-set:
	for img in ${TEST_IMAGES} ; do \
          ${MAKE} TESTARGS="${TESTARGS}" IMGARGS="$$img" clean-single ; \
        done
	rm -fv ${OT_STRIP_IMAGE}

######################################################################
# Submake for test on single image -- caller must define ${IMGARGS}

IMGARGS := 

ifneq "/${IMGARGS}" "/"

IMGFIELDS := ${subst :, ,${IMGARGS}}
IMG_NAME :=  ${word 1,${IMGFIELDS}}
IMG_EXT :=   ${word 2,${IMGFIELDS}}

IN_IMAGE := ${IN_DIR}/${IMG_NAME}.${IMG_EXT}

OT_IMAGE := ${OT_PREFIX}-${IMG_NAME}.${IMG_EXT}

test-single: ${OT_IMAGE}

${OT_IMAGE}: ${IN_IMAGE} ${PROGDIR}/${PROG} Makefile
	mkdir -p ${OT_DIR}/${FILTER_NAME}
	cat ${IN_IMAGE} \
          | time ${PROGDIR}/${PROG} \
              -weights ${WT_IMAGE} \
              -filter ${FILTER_SPECS} \
              -verbose \
          > $@
ifeq "/${SHOW}" "/YES"
	-if [ -s ${OT_IMAGE} ] ; then \
          display -title '%d/%f' ${OT_IMAGE} ${IN_IMAGE} ; \
        fi
endif  
ifeq "/${GPROF}" "/YES"
	-if [[ -r gmon.out ]] ; then \
          gprof \
            -I .. \
            -I ../../../JSLIBS/libimg \
            -I ../../../JSLIBS/libjs \
            --flat-profile \
            ${PROGDIR}/${PROG} \
            gmon.out ; \
        fi
endif  

clean-single:
	rm -fv ${OT_IMAGE}

endif
# End of ${IMGARGS}
######################################################################

endif
# End of ${TESTARGS}
######################################################################
