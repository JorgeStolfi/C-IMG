# Tests of the pnmprojmap tool
# Last edited on 2023-10-15 12:25:16 by stolfi

PROG := pnmprojmap
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

DEBUG := -verbose -debug 140 195

# Fields: {MAP}:{IMG}:{EXT}:{XAXIS}:{YAXIS}:{ORG}:{IUNIT}:{OUNIT}:{NXO}:{NYO}
# where {ORG} is "C" to use "center", "O" to use "org 20 30", "D" for the default,
# and "{NXO}:{NYO} is the output image size, or "DEF:DEF" for the default.

TESTS_AXES := \
  ident:ere-1:ppm:r:d:C:1.000:1.000:500:700 \
  ident:ere-1:ppm:r:u:C:1.000:1.000:500:700 \
  ident:ere-1:ppm:l:d:C:1.000:1.000:500:700 \
  ident:ere-1:ppm:l:u:C:1.000:1.000:500:700

TESTS_ORIGIN := \
  ident:ere-1:ppm:r:d:O:1.000:1.000:500:700 \
  ident:ere-1:ppm:r:u:O:1.000:1.000:500:700 \
  ident:ere-1:ppm:l:d:O:1.000:1.000:500:700 \
  ident:ere-1:ppm:l:u:O:1.000:1.000:500:700 \
   \
  ident:ere-1:ppm:r:d:D:1.000:1.000:500:700 \
  ident:ere-1:ppm:r:u:D:1.000:1.000:500:700 \
  ident:ere-1:ppm:l:d:D:1.000:1.000:500:700 \
  ident:ere-1:ppm:l:u:D:1.000:1.000:500:700

TESTS_ISCALE := \
  ident:ere-1:ppm:r:d:C:0.500:1.000:500:700 \
  ident:ere-1:ppm:r:u:C:0.500:1.000:500:700 \
  ident:ere-1:ppm:l:d:C:0.500:1.000:500:700 \
  ident:ere-1:ppm:r:d:C:0.500:1.000:500:700 \
   \
  ident:ere-1:ppm:r:d:O:0.500:1.000:500:700
 
TESTS_OSCALE := \
  ident:ere-1:ppm:r:d:C:1.000:3.000:500:700 \
  ident:ere-1:ppm:r:u:C:1.000:3.000:500:700 \
  ident:ere-1:ppm:l:d:C:1.000:3.000:500:700 \
  ident:ere-1:ppm:r:d:C:1.000:3.000:500:700 \
   \
  ident:ere-1:ppm:r:d:O:1.000:3.000:500:700
  
TESTS_DEFAULT_OSIZE := \
  ident:abi-1:ppm:r:d:C:1.000:1.000:DEF:DEF \
  ident:flo-1:ppm:r:d:O:0.500:3.000:DEF_DEF
 
TESTS_BASIC_MAPS := \
  rot90:ere-1:ppm:r:d:C:1.000:4.000:500:700 \
  shr45:ere-1:ppm:r:d:C:1.000:4.000:500:700 \
  rot45:ere-1:ppm:r:d:C:1.000:4.000:500:700 \
  shf05:ere-1:ppm:r:d:C:1.000:4.000:500:700

TESTS_GRAY := \
  ident:ere-1:pgm:r:d:C:1.000:1.000:500:700 \
  shf05:ere-1:pgm:r:d:C:1.000:1.000:500:700 \
  cubic:ere-1:pgm:r:d:C:1.000:1.000:500:700

TESTS := ${TESTS_AXES} ${TESTS_ORIGIN} ${TESTS_OSCALE} ${TESTS_ISCALE} ${TESTS_BASIC_MAPS} ${TESTS_GRAY}
# TESTS := ${TESTS_BASIC_MAPS} ${TESTS_OSCALE} ${TESTS_ISCALE}
# TESTS := ${TESTS_BASIC_MAPS}
# TESTS := ${TESTS_OSCALE} ${TESTS_ISCALE}
 
MAKEFILE := Makefile

.PHONY:: \
  run-tests single-test map-single compare-single

all: clean run-tests

run-tests:
	for tt in ${TESTS}; do \
          ${MAKE} SPECS=$$tt single-test; \
        done

clean::
	rm -f out/*.p?m

SHOW := YES

########################################################################
# For recursive "make" of a single test. Caller must define ${SPECS}.

ifneq "/${SPECS}" "/"

SPECFIELDS := ${subst :, ,${SPECS}}
MAP :=    ${word  1,${SPECFIELDS}}
IMG :=    ${word  2,${SPECFIELDS}}
EXT :=    ${word  3,${SPECFIELDS}}
XAXIS :=  ${word  4,${SPECFIELDS}}
YAXIS :=  ${word  5,${SPECFIELDS}}
ORG :=    ${word  6,${SPECFIELDS}}
IUNIT :=  ${word  7,${SPECFIELDS}}
OUNIT :=  ${word  8,${SPECFIELDS}}
NXO :=    ${word  9,${SPECFIELDS}}
NYO :=    ${word 10,${SPECFIELDS}}

UNIT_TAG := ${subst .,,${IUNIT}}-${subst .,,${OUNIT}}
USYS_TAG := ${XAXIS}${YAXIS}${ORG}-${UNIT_TAG}

# Input image file:
IMG_IN_FILE := data/${IMG}.${EXT}
MAP_IN_FILE := data/${MAP}.map

# Input image attributes:
IMG_ATTRS := 
ifeq "${ORG}" "C"
  IMG_ATTRS := ${IMG_ATTRS} center
else ifeq "${ORG}" "O"
  IMG_ATTRS := ${IMG_ATTRS} org 30 40
endif
ifneq "${IUNIT}" "1.000"
  IMG_ATTRS := ${IMG_ATTRS} unit ${IUNIT}
endif

# Assumed input size:
ISIZE_PIX := 350 490
ISIZE_USR := ${shell echo ${ISIZE_PIX} | gawk '//{ printf "%03d %03d\n", $$1/${IUNIT}, $$2/${IUNIT} }'}
ifeq "${NXO}:${NYO}" "DEF:DEF"
  # Default output size:
  OSIZE_OP := 
  OSIZE_USR := ${ISIZE_USR}
  OSIZE_TAG := DEFxDEF
else
  OSIZE_OP := -oSize ${NXO} ${NYO}
  OSIZE_USR := ${NXO} ${NYO}
  OSIZE_TAG := ${NXO}x${NYO}
endif
OSIZE_PIX := ${shell echo ${OSIZE_USR} | gawk '//{ printf "%03d %03d\n", $$1*${OUNIT}, $$2*${OUNIT} }'}

INTER_OP := -interpolate 1
INTER_TAG := int1

# Prefix for all output and reference file names of this test:
TEST_PREF := ${IMG}-${MAP}-${USYS_TAG}-${OSIZE_TAG}-${INTER_TAG}

# Output image and map files:
IMG_OT_FILE := out/${TEST_PREF}-${IMG}.${EXT}
MAP_OT_FILE := out/${TEST_PREF}-${MAP}.map

# Reference files for comparison with output:
IMG_EX_FILE := ${subst out/,ref/,${IMG_OT_FILE}}
MAP_EX_FILE := ${subst out/,ref/,${MAP_OT_FILE}}

IMG_OT_EX_DIF_FILE := out/${TEST_PREF}-dif.${EXT}
MAP_OT_EX_DIF_FILE := out/${TEST_PREF}-dif.txt

single-test: map-single compare-single

map-single: ${IMG_OT_FILE}
	touch out/.tried-mapping
        
compare-single: ${IMG_OT_EX_DIF_FILE} ${MAP_OT_EX_DIF_FILE}

${IMG_OT_EX_DIF_FILE}: out/.tried-mapping
	../compare_mapped_image.sh ${IMG_OT_FILE} ${IMG_EX_FILE} ${IMG_OT_EX_DIF_FILE} ${SHOW}

${MAP_OT_EX_DIF_FILE}: out/.tried-mapping
	../compare_output_points_map.sh ${MAP_OT_FILE} ${MAP_EX_FILE} ${MAP_OT_EX_DIF_FILE} ${SHOW}

${IMG_OT_FILE}: ${IMG_IN_FILE} ${MAP_IN_FILE} ${PROGDIR}/${PROG}
	@echo "ISIZE_PIX = ${ISIZE_PIX}" 1>&2
	@echo "ISIZE_USR = ${ISIZE_USR}" 1>&2
	@echo "OSIZE_USR = ${OSIZE_USR}" 1>&2
	@echo "OSIZE_PIX = ${OSIZE_PIX}" 1>&2
	${PROGDIR}/${PROG} \
            -inPrefix "data/" \
            -xAxis ${XAXIS} \
            -yAxis ${YAXIS} \
            -fromMatrix \
            -image ${IMG} ${EXT} ${IMG_ATTRS} matrix ${MAP} \
            -oCenter \
            -oUnit ${OUNIT} \
            ${OSIZE_OP} \
            ${INTER_OP} \
            -outPrefix "out/${TEST_PREF}-" \
            ${DEBUG}

endif
# End  ${SPECS} section.
########################################################################

