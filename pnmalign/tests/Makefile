# Tests of the pnmalign tool
# Last edited on 2011-06-05 23:51:52 by stolfi

PROG := pnmalign
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  gsm-pgm-3:10:10:gsm-1,gsm-2,gsm-3:pgm

NONTESTS := \
  abi-ppm-3:40:20:abi-1,abi-2,abi-3:ppm \
  grd-pgm-3:20:15:grd-1,grd-2,grd-3:pgm \
  grd-ppm-3:20:15:grd-1,grd-2,grd-3:ppm

.PHONY:: \
  single-test compare-images \
  ident rot45 persp 

all:
	for tt in ${TESTS}; do \
          ${MAKE} SPECS=$$tt single-test; \
        done

clean:
	rm -f core out/*.p?m *-ref.p?m *-dif.p?m 

ifneq "/${SPECS}" "/"

COMMA := ,
SPECFIELDS := ${subst :, ,${SPECS}}
TST := ${word 1,${SPECFIELDS}}
DMX := ${word 2,${SPECFIELDS}}
RAD := ${word 3,${SPECFIELDS}}
LST := ${subst ${COMMA}, ,${word 4,${SPECFIELDS}}}
EXT := ${word 5,${SPECFIELDS}}

INDISP := ${TST}-inp.dsp
INFILES := ${addsuffix .${EXT},${LST}}

OTDISP := out/${TST}.dsp
OTFILES := ${addsuffix -shf.${EXT},${LST}}

single-test: ${OTDISP}
	cat ${OTDISP}
	transform-images ${OTDISP} ${INFILES}

INARGS := ${shell make-align-args ${INDISP} ${INFILES}}

${OTDISP}: Makefile ${INFILES} ${PROGDIR}/${PROG}
	${PROGDIR}/${PROG} -radius ${RAD} -dmax ${DMX} \
            ${INARGS} \
          > ${OTDISP}

endif
