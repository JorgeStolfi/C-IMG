# Tests of the pnmprojmap tool
# Last edited on 2023-10-15 12:43:29 by stolfi

PROG := pnmprojmap
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DEBUG := -verbose -debug 140 195

# Fields: {IMG_PR}:{XAXIS}:{YAXIS}:{IUNIT}:{OUNIT}
# where {IMG_PR} is the image name minus the "-{DIGIT}" part.

TESTS_KMNO4 := \
  kmn:r:d:1.000:0.500 \
  kmn:r:d:0.500:3:000

TESTS := ${TESTS_KMNO4}

MAKEFILE := Makefile

.PHONY:: \
  run-tests single-test map-images compare-single

all: clean run-tests

run-tests:
	for tt in ${TESTS}; do \
          ${MAKE} SPECS=$$tt single-test; \
        done

clean::
	rm -f out/*.p?m

SHOWIMG := display

########################################################################
# For recursive "make" of a single test. Caller must define ${SPECS}.

ifneq "/${SPECS}" "/"

# {IMG_PR}:{XAXIS}:{YAXIS}:{ORG}:{IUNIT}:{OUNIT}:{NXO}:{NYO}

SPECFIELDS := ${subst :, ,${SPECS}}
IMG_PR := ${word  1,${SPECFIELDS}}
XAXIS :=  ${word  2,${SPECFIELDS}}
YAXIS :=  ${word  3,${SPECFIELDS}}
IUNIT :=  ${word  4,${SPECFIELDS}}
OUNIT :=  ${word  5,${SPECFIELDS}}

UNIT_TAG := ${subst .,,${IUNIT}}-${subst .,,${OUNIT}}
USYS_TAG := ${XAXIS}${YAXIS}-${UNIT_TAG}

# The first one will be the reference image:
EXT := ppm

# Input files:
IMG_IN_FILES := ${wildcard data/${IMG_PR}-[0-9].${EXT}}
PTS_IN_FILES := ${add_suffix .pts,${basename IMG_IN_FILES}}
IMG_NAMES := ${basename ${notdir ${IMG_IN_FILES}}}

# Input image attributes:
IMG_ATTRS := :${EXT}
ifeq "${ORG}" "C"
  IMG_ATTRS := ${IMG_ATTRS}:center
else ifeq "${ORG}" "O"
  IMG_ATTRS := ${IMG_ATTRS}:org:30:40
endif
ifneq "${IUNIT}" "1.000"
  IMG_ATTRS := ${IMG_ATTRS}:unit:${IUNIT}
endif

IMAGE_OPS := ${subst :, ,${addprefix -image:,${addsuffix ${IMG_ATTRS},${IMG_NAMES}}}}

# Output image size:
NXO := 500
NYO := 700
OSIZE_TAG := ${NXO}x${NYO}

# Prefix for all output and reference file names of this test:
TEST_PREF := ${IMG_PR}-${USYS_TAG}-${OSIZE_TAG}

single-test: map-images
	for img in ${IMG_NAMES} do \
          ${MAKE} IMG_NAME=$$img compare-single ; \
        done
        
map-images: ${IMG_IN_FILES} ${PTS_IN_FILES} ${PROGDIR}/${PROG}
	@ls -l data/${IMG_PR}-[0-9].${EXT} 1>&2
	@echo "IMG_NAMES = ${IMG_NAMES}" 1>&2
	@echo "IMG_IN_FILES = ${IMG_IN_FILES}" 1>&2
	mkdir -v -p out
	${PROGDIR}/${PROG} \
            -inPrefix "data/" \
            -xAxis ${XAXIS} \
            -yAxis ${YAXIS} \
            -fromPoints \
            ${IMAGE_OPS} \
            -oCenter \
            -oUnit ${OUNIT} \
            -oSize ${NXO} ${NYO} \
            -interpolate 1 \
            -outPrefix "out/${TEST_PREF}-"
	out/.tried-mapping
        
######################################################################
# For recursive "make" of a single image in a single test. Caller must define ${IMG_NAME}.
ifneq "/IMG_NAME" "/"

# Input image and point files
IMG_IN_FILE := data/${IMG_NAME}.${EXT}
PTS_IN_FILE := ${add_suffix .pts,${basename IMG_IN_FILE}}
  
# Output image, points, and map files:
IMG_OT_FILE := ${subst data/,out/${TEST_PREF}-,${IMG_IN_FILE}}
PTS_OT_FILE := ${subst .${EXT},.pts,${IMG_OT_FILE}}
MAP_OT_FILE := ${subst .${EXT},.map,${IMG_OT_FILE}}
 
compare-single: ${IMG_OT_EX_DIF_FILE} ${PTS_OT_EX_DIF_FILE} ${MAP_OT_EX_DIF_FILE}

${IMG_OT_EX_DIF_FILE}: out/.tried-mapping
	../compare_mapped_image.sh ${IMG_OT_FILE} ${IMG_EX_FILE} ${IMG_OT_EX_DIF_FILE} ${SHOW}

${PTS_OT_EX_DIF_FILE}: out/.tried-mapping
	../compare_output_points_map.sh ${PTS_OT_FILE} ${PTS_EX_FILE} ${PTS_OT_EX_DIF_FILE} ${SHOW}

${MAP_OT_EX_DIF_FILE}: out/.tried-mapping
	../compare_output_points_map.sh ${MAP_OT_FILE} ${MAP_EX_FILE} ${MAP_OT_EX_DIF_FILE} ${SHOW}
        
endif
# End  ${IMG_NAME} section.
######################################################################

endif
# End  ${SPECS} section.
########################################################################

