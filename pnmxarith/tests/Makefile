# Last edited on 2009-01-06 17:52:31 by stolfi

PROG := pnmxarith
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  ramp:ppm:-add \
  ramp:ppm:-add%-scale%0.5 \
  ramp:ppm:-multiply
  
all:
	for t in ${TESTS}; do \
          ${MAKE} TEST=$$t clean-single single ; \
        done
        
clean::
	for t in ${TESTS} ${NONTESTS}; do \
          ${MAKE} TEST=$$t clean-single ; \
        done

SHOW := display -title '%f'

######################################################################
# Sub-make - caller must define ${TEST}

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

FIELDS := ${subst :, ,${TEST}}
NAME := ${word 1, ${FIELDS}}
EXT :=  ${word 2, ${FIELDS}}
OPS :=  ${subst %, ,${word 3, ${FIELDS}}}

single: ${NAME}-r.${EXT} 
	${SHOW} -filter Point -resize '200%' \
          ${NAME}-0.${EXT} \
          ${NAME}-1.${EXT} \
          ${NAME}-r.${EXT}

${NAME}-r.${EXT}: ${PROGDIR}/${PROG} ${NAME}-0.${EXT} ${NAME}-1.${EXT}
	${PROGDIR}/${PROG} \
            ${OPS} \
            ${NAME}-0.${EXT} \
            ${NAME}-1.${EXT} \
          > ${NAME}-r.${EXT}
	identify ${NAME}-r.${EXT}

clean-single:
	-rm -f ${NAME}-r.${EXT}

endif
# end ${TEST} sub-make.
######################################################################

######################################################################
# For the "ramp" test:

MAKERAMP := ppmxramp

ifeq "/${EXT}" "/ppm" 
  BLACK := rgbi:0.0/0.0/0.0
  BEIGE := rgbi:1.0/1.0/0.5
  LILAC := rgbi:1.0/0.5/1.0
endif
ifeq "/${EXT}" "/pgm" 
  BLACK := 0.0
  BEIGE := 0.4
  LILAC := 0.7
endif

ramp-0.${EXT}: ${MAKERAMP}
	${MAKERAMP} \
            256 256 \
            000 000 ${BLACK} \
            255 000 ${BEIGE} \
            255 255 ${BEIGE} \
          > ramp-0.${EXT}
	${SHOW} ramp-0.${EXT}

ramp-1.${EXT}: ${MAKERAMP}
	${MAKERAMP} \
            256 256 \
            000 000 ${BLACK} \
            000 255 ${LILAC} \
            255 255 ${LILAC} \
          > ramp-1.${EXT}
	${SHOW} ramp-1.${EXT}
          
# end creation of images for "ramp" test.
######################################################################
