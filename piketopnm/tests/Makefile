# Last edited on 2023-02-04 07:07:28 by stolfi

PROG := piketopnm

PROGDIR  := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# ${IMGNAME}:${DEB_PAD}:${DEB_INT}:${DEB_SQZ}:${EXTR_CHN}:${RC}:${GC}:${BC}:${DEN}

TESTS := \
  led-2:Y:N:N:--:20966:13983:18984:12000 \
  led-2:N:Y:N:--:20966:13983:18984:12000 \
  led-2:N:N:Y:--:20966:13983:18984:12000 \

NON_TESTS := \
  syn-1:N:N:N:--:10000:10000:10000:10000 \
  syn-1:N:N:N:G0:10000:10000:10000:10000 \
  syn-1:N:N:N:R0:10000:10000:10000:10000 \
  syn-1:N:N:N:B0:10000:10000:10000:10000 \
  syn-1:N:N:N:G1:10000:10000:10000:10000 \
  \
  led-1:N:N:N:--:10000:10000:10000:10000 \
  led-1:Y:N:N:--:10000:10000:10000:10000 \
  led-1:Y:N:N:--:20966:13983:18984:12000 \
  led-1:N:Y:N:--:20966:13983:18984:12000 \
  hal-1:N:Y:N:--:10000:10000:10000:10000 \
  \
  led-1:N:N:N:G0:10000:10000:10000:10000 \
  led-1:N:N:N:R0:10000:10000:10000:10000 \
  led-1:N:N:N:B0:10000:10000:10000:10000 \
  led-1:N:N:N:G1:10000:10000:10000:10000 \
  \
  led-2:N:N:N:--:10000:10000:10000:10000 \
  led-2:N:N:N:G0:10000:10000:10000:10000 \
  led-2:N:N:N:R0:10000:10000:10000:10000 \
  led-2:N:N:N:B0:10000:10000:10000:10000 \
  led-2:N:N:N:G1:10000:10000:10000:10000 \
  
.PHONY:: run-all run-single

all: run-all

clean:: 
	( cd out && rm -f *.ppm *.pgm )

run-all: ; \
  for tst in ${TESTS}; do \
    ${MAKE} TEST_ARGS=$${tst} run-single; \
  done

######################################################################
# Recursive section -- client must define ${TEST}
TEST_ARGS := 
ifneq "/${TEST_ARGS}" "/"

TEST_FIELDS := ${subst :, ,${TEST_ARGS}}
IMGNAME   := ${word 1, ${TEST_FIELDS}}
DEB_PAD   := ${word 2, ${TEST_FIELDS}}
DEB_INT   := ${word 3, ${TEST_FIELDS}}
DEB_SQZ   := ${word 4, ${TEST_FIELDS}}
EXTR_CHN  := ${word 5, ${TEST_FIELDS}}
RC        := ${word 6, ${TEST_FIELDS}}
GC        := ${word 7, ${TEST_FIELDS}}
BC        := ${word 8, ${TEST_FIELDS}}
DEN       := ${word 9, ${TEST_FIELDS}}

EXT := EXT.IS.UNDEFINED

DEBAYER_OPS := 
ifeq "/${DEB_PAD}" "/Y"
  EXT := ppm
  DEBAYER_OPS := -debayer pad
endif
ifeq "/${DEB_INT}" "/Y"
  EXT := ppm
  DEBAYER_OPS := -debayer interpolate
endif
ifeq "/${DEB_SQZ}" "/Y"
  EXT := pgm
  DEBAYER_OPS := -debayer squeeze
endif
ifneq "/${EXTR_CHN}" "/--"
  EXT := pgm
  DEBAYER_OPS := -extract ${subst -,,${EXTR_CHN}}
endif

SCRIPT := data/${IMGNAME}.doit
INFILE := data/${IMGNAME}.raw
OTFILE := out/${IMGNAME}-${DEB_PAD}-${DEB_INT}-${DEB_SQZ}-${EXTR_CHN}-${RC}-${GC}-${BC}-${DEN}.${EXT}

SHOWIMG := display
PSVIEW := evince
# PSVIEW := gv

run-single: ${OTFILE}

${OTFILE}: ${PROGDIR}/${PROG} ${INFILE}
	@echo INFILE = ${INFILE}
	${PROGDIR}/${PROG} \
            ${DEBAYER_OPS} \
            -balance ${RC} ${GC} ${BC} / ${DEN} \
            -verbose \
          < ${INFILE} \
          | pnmnoraw \
          > ${OTFILE}
	-${SHOWIMG} ${OTFILE}

endif
# End ${TEST_ARGS} section
######################################################################
