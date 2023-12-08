# Tests of the pnmgtran tool
# Last edited on 2010-07-02 16:56:28 by stolfi

PROG := pnmfftfilter
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# W_MIN_X:W_MIN_Y:W_MAX_X:W_MAX_Y:ORIG:EXT

SINGLE_TEST := \
  -pass:0020.000:0020.000:0100.000:0100.000:ula-1:ppm

RND_TESTS := \
  -pass:0000.000:0000.000:0002.500:0002.500:rnd-1:pgm \
  -pass:0002.500:0002.500:0005.000:0005.000:rnd-1:pgm \
  -pass:0005.000:0005.000:0010.000:0010.000:rnd-1:pgm \
  -pass:0010.000:0010.000:0020.000:0020.000:rnd-1:pgm \
  -pass:0020.000:0020.000:0040.000:0040.000:rnd-1:pgm \
  -pass:0040.000:0040.000:0080.000:0080.000:rnd-1:pgm \
  -pass:0080.000:0080.000:0160.000:0160.000:rnd-1:pgm \
  -pass:0160.000:0160.000:0320.000:0320.000:rnd-1:pgm \
  -pass:0320.000:0320.000:0640.000:0640.000:rnd-1:pgm \
  -pass:0640.000:0640.000:9999.000:9999.000:rnd-1:pgm \
  -pass:9999.000:9999.000:infinity:infinity:rnd-1:pgm \
  -pass:0010.000:0100.000:0030.000:0300.000:rnd-1:pgm \
  -kill:0100.000:0100.000:0300.000:0300.000:rnd-1:pgm
  
ULA_W_ZER := 0000.000
ULA_W_MIN := 0008.000
ULA_W_MAX := 0050.000
ULA_W_INF := infinity

ULA_TESTS := \
  -pass:${ULA_W_ZER}:${ULA_W_ZER}:${ULA_W_MIN}:${ULA_W_MIN}:ula-m:ppm \
  -pass:${ULA_W_MIN}:${ULA_W_MIN}:${ULA_W_MAX}:${ULA_W_MAX}:ula-m:ppm \
  -kill:${ULA_W_MIN}:${ULA_W_MIN}:${ULA_W_MAX}:${ULA_W_MAX}:ula-m:ppm \
  -pass:${ULA_W_MAX}:${ULA_W_MAX}:${ULA_W_INF}:${ULA_W_INF}:ula-m:ppm

NONTESTS_GRAY := \
  -pass:0020.000:0020.000:0060.000:0060.000:flo-1:pgm

NONTESTS_COLOR := \
  -pass:0020.000:0020.000:0060.000:0060.000:flo-1:ppm \
  \
  -pass:0000.000:0000.000:9999.000:9999.000:abi-1:ppm \
  -pass:0000.000:0000.000:0020.000:0020.000:abi-1:ppm \
  -pass:0020.000:0020.000:0060.000:0060.000:abi-1:ppm \
  -pass:0060.000:0060.000:9999.000:9999.000:abi-1:ppm \
  \
  -pass:0000.000:0000.000:0020.000:0020.000:ula-1:ppm \
  -pass:0020.000:0020.000:0100.000:0100.000:ula-1:ppm \
  -pass:0100.000:0100.000:9999.000:9999.000:ula-1:ppm

# TESTS := ${SINGLE_TEST}
# TESTS := ${RND_TESTS}
TESTS := ${ULA_TESTS}

MAKEFILE := Makefile

.PHONY:: \
  run-tests single-test compare-images

all: clean run-tests

run-tests:
	for tt in ${TESTS}; do \
          ${MAKE} SPECS=$$tt single-test; \
        done

clean::
	rm -f *-out.p?m *-ref.p?m *-dif.p?m 

SHOWIMG := display -title '%f'
# SHOWIMG := xv

########################################################################
# For recursive "make" of a single test. Caller must define ${SPECS}.

ifneq "/${SPECS}" "/"

SPECFIELDS := ${subst :, ,${SPECS}}
FKIND :=    ${word 1,${SPECFIELDS}}
W_MIN_X :=  ${word 2,${SPECFIELDS}}
W_MIN_Y :=  ${word 3,${SPECFIELDS}}
W_MAX_X :=  ${word 4,${SPECFIELDS}}
W_MAX_Y :=  ${word 5,${SPECFIELDS}}
IMG :=    ${word 6,${SPECFIELDS}}
EXT :=    ${word 7,${SPECFIELDS}}

DEBUG := -verbose

INFILE := data/${IMG}.${EXT}
OTFILE := out/${IMG}-${FKIND}-${W_MIN_X}-${W_MIN_Y}--${W_MAX_X}-${W_MAX_Y}-out.${EXT}

${OTFILE}: ${INFILE} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} \
            ${FKIND} \
            -from ${W_MIN_X} ${W_MIN_Y} \
            -to ${subst infinity,+oo,${W_MAX_X}} ${subst infinity,+oo,${W_MAX_Y}} \
            ${DEBUG} \
            -range AUTO \
            -maxval 10000 \
            ${INFILE} \
          > ${OTFILE}

single-test: ${OTFILE} ${INFILE}
	${SHOWIMG} ${OTFILE} ${INFILE}

endif
# End ${SPECS}.
########################################################################
