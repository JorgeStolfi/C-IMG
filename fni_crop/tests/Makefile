# Last edited on 2024-12-25 07:17:32 by stolfi

SHELL =	/bin/bash

PROG := fni_crop
        
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# {LEFT}:{RIGHT}:{WIDTH}:{TOP}:{BOTTOM}:{HEIGHT}:{YAXIS}

ALL_TESTS := \
  NONE:NONE:NONE:NONE:NONE:NONE:d \
  NONE:NONE:NONE:NONE:NONE:NONE:u \
   \
  +150:NONE:NONE:+220:NONE:NONE:d \
  NONE:+150:NONE:NONE:+220:NONE:d \
  NONE:NONE:0150:NONE:NONE:0220:d \
   \
  +150:+010:NONE:+220:+040:NONE:d \
  NONE:+150:0150:NONE:+120:0120:d \
  +150:NONE:0150:+120:NONE:0120:d \
   \
  -050:NONE:NONE:-020:NONE:NONE:d \
  NONE:-050:NONE:NONE:-020:NONE:d \
   \
  +150:+150:0200:+120:+050:0180:d \
  +150:+150:0200:+120:+050:0180:u \

TESTS := ${ALL_TESTS}

.PHONY:: runtest-all runtest

all: runtest-all

# ====================================================================
# Cleanup

clean::
	-/bin/rm -f out/*.fni

# ====================================================================
# Running tests

runtest-all: 
	-mkdir -p out
	-( cd out && rm -f *.fni *.ppm )
	for tt in ${TESTS} ; do \
          ${MAKE} TEST=$$tt runtest ; \
        done

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"
#######################################################################
# Recursive "make" to create the test ${TEST}

ARGS := ${subst :, ,${TEST}}
LEFT    := ${word 1,${ARGS}}
RIGHT   := ${word 2,${ARGS}}
WIDTH   := ${word 3,${ARGS}}
TOP     := ${word 4,${ARGS}}
BOTTOM  := ${word 5,${ARGS}}
HEIGHT  := ${word 6,${ARGS}}
YAXIS   := ${word 7,${ARGS}}

OPTS := -yAxis ${YAXIS}

ifneq "${LEFT}" "NONE"
OPTS += -left ${LEFT}
endif

ifneq "${RIGHT}" "NONE"
OPTS += -right ${RIGHT}
endif

ifneq "${WIDTH}" "NONE"
OPTS += -width ${WIDTH}
endif

ifneq "${TOP}" "NONE"
OPTS += -top ${TOP}
endif

ifneq "${BOTTOM}" "NONE"
OPTS += -bottom ${BOTTOM}
endif

ifneq "${HEIGHT}" "NONE"
OPTS += -height ${HEIGHT}
endif

INFILE_FNI := in/frutas01.fni
OTFILE_FNI := out/frutasgd:${LEFT}:${RIGHT}:${WIDTH}:${TOP}:${BOTTOM}:${HEIGHT}:${YAXIS}.fni
OTFILE_FNI := ${subst -,m,${OTFILE_FNI}}
OTFILE_FNI := ${subst +,p,${OTFILE_FNI}}
OTFILE_FNI := ${subst :,-,${OTFILE_FNI}}
OTFILE_PPM := ${subst .fni,.ppm,${OTFILE_FNI}}

runtest: ${PROGDIR}/${PROG} ${INFILE_FNI} Makefile
	@echo "OPTS = ${OPTS}"
	${PROGDIR}/${PROG} \
          ${OPTS} \
          -verbose \
          -in ${INFILE_FNI} \
          -out ${OTFILE_FNI}
	fni_to_pnm -yAxis down < ${OTFILE_FNI} > ${OTFILE_PPM}
	( display -title '%f' -filter Box -resize '200%' ${OTFILE_PPM} )

# 
# End make ${ZFUNC}, ${SIGMA}
######################################################################
endif
