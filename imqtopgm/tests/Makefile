# Last edited on 2023-02-07 22:54:34 by stolfi

PROG := imqtopgm
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

IMAGES := \
  F429A13 \
  F427A03 \
  C2064910 \
  C2065022
  
all:
	mkdir -p out
	rm -f out/*pgm
	for img in ${IMAGES} ; do \
          ${MAKE} IMG=$${img} single ; \
        done
        
# ######################################################################
# Recursive section, user must define ${IMG}
ifneq "/${IMG}" "/"

IN_FILE := data/${IMG}.IMQ
OUT_FILE := out/${IMG}.pgm
REF_FILE := out/${IMG}-ref.pgm

single: create display

create: ${OUT_FILE}

${OUT_FILE}: ${PROGDIR}/${PROG} ${IN_FILE}
	../orig-compile/vdcomp < ${IN_FILE} > ${REF_FILE}
	${PROGDIR}/${PROG} < ${IN_FILE} > ${OUT_FILE}
        
display: ${OUT_FILE}
	display -title '%f' ${OUT_FILE}

endif
# end ${IMG} section
# ######################################################################
        
clean::

