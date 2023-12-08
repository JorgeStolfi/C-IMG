# Last edited on 2023-09-24 16:08:37 by stolfi

PROG := pnmdiffop
PROGDIR := ..

SHELL := /bin/bash

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

SHOW_SINGLE := NO
SHOW_STRIP := YES
GPROF := NO

# Format: "{IMG_NAME}:{FORMAT}"

TEST_IMAGES_SOME_GRAY := \
  bullseye-128:pgm \
  cuneif01-128:pgm \
  branches-128:pgm \
  nerays01-128:pgm \
  seasqirt-128:pgm \
  frutas02-128:pgm

TEST_IMAGES_SOME_COLOR := \
  frutas01-128:ppm
 
TEST_IMAGES_ALL_GRAY := \
  doctyp01-128:pgm \
  bullseye-128:pgm \
  nerays01-128:pgm \
  ballmou1-128:pgm \
  branches-128:pgm \
  bricks02-128:pgm \
  cuneif01-128:pgm \
  seasqirt-128:pgm \
  frutas02-128:pgm
 
TEST_IMAGES_ALL_COLOR := \
  charsvms-128:ppm \
  frutas01-128:ppm \
  maiden01-128:ppm
 
# TEST_IMAGES := ${TEST_IMAGES_ALL_GRAY}
# TEST_IMAGES := ${TEST_IMAGES_SOME_GRAY}
TEST_IMAGES := ${TEST_IMAGES_ALL_COLOR}

# Format: "{SMOOTHABLE}:{OPNAME}"

SOME_OPERATORS := \
  0:variance \
  0:deviation \
  1:laplacian \

ALL_OPERATORS := \
  0:average \
  0:variance \
  0:deviation \
  0:f \
  1:fx \
  1:fy \
  1:fxx \
  1:fxy \
  1:fyy \
  0:fxxy \
  0:fxyy \
  0:fxxyy \
  1:gradient \
  1:laplacian \
  0:orthicity \
  0:elongation \
  0:linf \
  0:linfx \
  0:linfy \
  0:linvar \
  0:lindev \

# OPERATORS := ${ALL_OPERATORS}
OPERATORS := ${SOME_OPERATORS}
#OPERATORS := 0:f
#OPERATORS := 1:fx
#OPERATORS := 1:fy
#OPERATORS := 1:gradient
#OPERATORS := 1:fxx

.PHONY:: \
  test-operator \
  test-set clean-set \
  test-single clean-single

# SQR_OPTION := Y
SQR_OPTION := N


ACTIONS := test-set
# ACTIONS := clean-set test-set

all:
	for op in ${OPERATORS}; do \
          ${MAKE} OPERATOR="$$op" test-operator ; \
        done

IN_DIR := in
OT_DIR := out

ifneq "/${OPERATOR}" "/"
######################################################################
# Submake for test of one operator on various images -- caller must define ${OPERATOR}

OPWDS := ${subst :, ,${OPERATOR}}
SMOOTHABLE := ${word 1, ${OPWDS}}
OPNAME :=     ${word 2, ${OPWDS}}

ifeq "/${SMOOTHABLE}" "/1"
  SMTS := N Y
else
  SMTS := N
endif

test-operator:
	for smt in ${SMTS} ; do \
          ${MAKE} \
            OPNAME="${OPNAME}" \
            SMT_OPTION="$$smt" \
            test-set ; \
        done

ifneq "/${SMT_OPTION}" "/"
######################################################################
# Submake for test of one operator with specific smoothed and squared options
# on various images -- caller must define ${SMT_OPTION} ${SQR_OPTION}

SMOOTHED_TAG := ${subst Y,smt,${subst N,raw,${SMT_OPTION}}}
SQUARED_TAG :=  ${subst Y,sqr,${subst N,uns,${SQR_OPTION}}}
SET_TAG := ${OPNAME}-${SMOOTHED_TAG}-${SQUARED_TAG}

IN_PARMS := ${IN_DIR}/${SET_TAG}.parms

PARMS := ${shell if [[ -s ${IN_PARMS} ]]; then cat ${IN_PARMS}; else echo "no ${IN_PARMS}, okay" 1>&2; fi}
PARMS_TAG := ${subst :,,${subst @,-,${subst .,,${subst +,p,${subst -,m,${PARMS}}}}}}

OT_PREFIX := ${OT_DIR}/${SET_TAG}/ot${PARMS_TAG}
OT_STRIP_PREFIX := ${OT_PREFIX}-strip
OT_STRIP_IMAGE := ${OT_STRIP_PREFIX}.png

OT_IMAGES := ${addprefix ${OT_PREFIX}-,${subst :,.,${TEST_IMAGES}}}

test-set:
	@echo "=== ${SET_TAG} ============================================================" 1>&2
	@echo "PARMS_TAG = ${PARMS_TAG}" 1>&2
	rm -fv ${OT_PREFIX}-*.p?m ${OT_PREFIX}-*.png
	rm -fv ${OT_STRIP_IMAGE}
	for img in ${TEST_IMAGES} ; do \
          ${MAKE} \
            OPNAME="${OPNAME}" \
            SMT_OPTION="${SMT_OPTION}" \
            SQR_OPTION="${SQR_OPTION}" \
            IMGARGS="$$img" test-single ; \
        done
	-@if [[ -s ${word 1, ${OT_IMAGES}} ]]; then \
  	  ls -l ${OT_IMAGES} ; \
          convert ${OT_IMAGES} \
            -bordercolor red -border 1x1 \
            +append ${OT_STRIP_IMAGE}; \
          if [[ ( "/${SHOW_STRIP}" == "/YES" ) && ( -s ${OT_STRIP_IMAGE} ) ]] ; then \
            display -title '%d/%f' -gamma 1.0 -filter Box -resize '200%' ${OT_STRIP_IMAGE} ; \
          fi; \
        else \
          echo "** no images generated!" 1>&2 ; \
	fi

clean-set:
	for img in ${TEST_IMAGES} ; do \
          ${MAKE} \
            OPNAME="${OPNAME}" \
            SMT_OPTION="${SMT_OPTION}" \
            SQR_OPTION="${SQR_OPTION}" \
            IMGARGS="$$img" clean-single ; \
        done
	rm -fv ${OT_STRIP_IMAGE}

ifneq "/${IMGARGS}" "/"
######################################################################
# Submake for test on single image -- caller must define ${IMGARGS}

IMGFIELDS :=    ${subst :, ,${IMGARGS}}
IMG_NAME :=     ${word 1,${IMGFIELDS}}
IMG_EXT :=      ${word 2,${IMGFIELDS}}

IN_IMAGE := ${IN_DIR}/${IMG_NAME}.${IMG_EXT}

OT_IMAGE := ${OT_PREFIX}-${IMG_NAME}.${IMG_EXT}

test-single: ${OT_IMAGE}

${OT_IMAGE}: ${IN_IMAGE} ${PROGDIR}/${PROG} Makefile
	@echo "... ${IMG_NAME} ......................................................" 1>&2
	@echo "IMGARGS = ${IMGARGS}" 1>&2
	@echo "OT_IMAGE = ${OT_IMAGE}" 1>&2
	mkdir -p ${OT_DIR}/${SET_TAG}
	cat ${IN_IMAGE} \
          | time ${PROGDIR}/${PROG} \
              -op ${OPNAME} \
              ${subst Y,-smoothed,${subst N,,${SMT_OPTION}}} \
              ${subst Y,-squared,${subst N,,${SQR_OPTION}}} \
              ${subst :, ,${subst @,-,${PARMS}}} \
              -verbose \
          > $@

ifeq "/${SHOW_SINGLE}" "/YES"
	-if [ -s ${OT_IMAGE} ] ; then \
          display -title '%d/%f' -gamma 1.0 ${OT_IMAGE} ${IN_IMAGE} ; \
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

# End of ${IMGARGS} 
######################################################################
endif

# End of ${SMT_OPTION}
######################################################################
endif

# End of ${OPERATOR}
######################################################################
endif
