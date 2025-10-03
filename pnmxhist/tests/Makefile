# Last edited on 2025-08-02 12:16:55 by stolfi

PROG := pnmxhist
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

PLOT := YES
SHOW := YES

# {IN_IMG}:{NEW_MAXVAL}:{CHANNEL}:{WT_IMG}:{IGNORE}
 
TESTS := \
  nypdfull.ppm:00255:B:NOWEIGHT.pgm:: \
  nypdfull.ppm:00255:B:nypdmask.pbm:: \
  nypdfull.ppm:00255:B:nypdmask.pgm::

TESTS_WITHOUT_MASK := \
  astolfi4.pbm:00001:Y:NOWEIGHT.pgm:: \
  astolfi4.pbm:00001:Y:NOWEIGHT.pgm:0: \
  bullseye.pgm:00128:Y:NOWEIGHT.pgm:: \
  bullseye.pgm:00032:Y:NOWEIGHT.pgm:1/3/5: \
  earthvis.pgm:65535:Y:NOWEIGHT.pbm:: \
  nypdfull.ppm:00255:B:NOWEIGHT.pgm:: \
  nypdfull.ppm:00255:R:NOWEIGHT.pgm::
  
TESTS_WITH_MASK := \
  astolfi4.pbm:00001:Y:astolfi4.pbm:: \
  bullseye.pgm:00128:Y:bullseye.pgm:: \
  earthvis.pgm:65535:Y:earthmsk.pbm:0: \
  earthvis.pgm:00255:Y:earthmsk.pbm:0: \
  nypdfull.ppm:00255:B:nypdmask.pbm:: \
  nypdfull.ppm:00255:B:nypdmask.pgm::
 
ALL_TESTS := \
  ${TESTS_WITHOUT_MASK} \
  ${TESTS_WITH_MASK}

.PHONY:: \
  test-single clean-single

all:
	for t in ${TESTS}; do \
          ${MAKE} TESTARGS="$$t" test-single ; \
        done

clean::
	for t in ${TESTS}; do \
          ${MAKE} TESTARGS="$$t" clean-single ; \
        done

######################################################################
# Submake for single test -- caller must define ${TESTARGS}

TESTARGS := 

ifneq "/${TESTARGS}" "/"

TESTFIELDS := ${subst :, ,${TESTARGS}}
IN_IMG :=     ${word 1,${TESTFIELDS}}
NEW_MAXVAL := ${word 2,${TESTFIELDS}}
CHANNEL :=    ${word 3,${TESTFIELDS}}
WT_IMG :=     ${word 4,${TESTFIELDS}}
IGNORE :=     ${word 5,${TESTFIELDS}}

EXT :=        ${extension ${IN_IMG}}

IN_DIR := in
WT_DIR := wt
OUT_DIR := out

# Input files:
IN_FILE := ${IN_DIR}/${IN_IMG}
WT_FILE := ${WT_DIR}/${WT_IMG}

# Program options:
ifeq "${suffix ${IN_IMG}}" ".ppm"
  CHANNEL_OPT := -channel ${CHANNEL}
else
  CHANNEL_OPT := 
endif

ifneq "${basename ${WT_IMG}}" "NOWEIGHT"
  MASK_OPT := -mask ${WT_FILE}
  DEP_FILES := ${WT_FILE} ${IN_FILE}
else
  MASK_OPT := 
  DEP_FILES := ${IN_FILE}
endif

ifneq "${IGNORE}" ""
  IGNORE_OPT := -ignore ${subst /, -ignore ,${IGNORE}}
else
  IGNORE_OPT := 
endif

# Tag for output file names:
IN_TAG := ${subst .,-,${IN_IMG}}
MV_TAG := m${NEW_MAXVAL}
WT_TAG := ${subst .,-,${WT_IMG}}
CH_TAG := c${CHANNEL}
IG_TAG := x${subst /,+,${IGNORE}}

OUT_TAG := ${IN_TAG}-${MV_TAG}-${CH_TAG}-${WT_TAG}-${IG_TAG}

OUT_PREFIX := ${OUT_DIR}/${OUT_TAG}

HIST_FILE := ${OUT_PREFIX}-hist.dat
PLOT_FILE := ${OUT_PREFIX}-hist.png

PLOT_SCRIPT := plot_pnm_hist.sh

ifeq "${PLOT}" "YES"
  TARGET_FILES := ${HIST_FILE} ${PLOT_FILE}
else
  TARGET_FILES := ${HIST_FILE} 
endif

test-single: clean-single ${TARGET_FILES}

${HIST_FILE}: ${PROGDIR}/${PROG} ${DEP_FILES}
	cat ${IN_FILE} \
        | pnmdepth ${NEW_MAXVAL} \
	| ${PROGDIR}/${PROG} \
            ${CHANNEL_OPT} \
            ${MASK_OPT} \
            ${IGNORE_OPT} \
            -verbose \
          > $@

${PLOT_FILE}: ${PLOT_SCRIPT} ${HIST_FILE}
	${PLOT_SCRIPT} ${HIST_FILE} ${PLOT_FILE} ${NEW_MAXVAL}
ifeq "/${SHOW}" "/YES"
	-if [ -s ${PLOT_FILE} ] ; then \
          display -title '%f' ${PLOT_FILE} ; \
        fi
endif  

clean-single:
	rm -f ${OUT_PREFIX}-*.*

endif
# End of ${TESTARGS}
######################################################################
