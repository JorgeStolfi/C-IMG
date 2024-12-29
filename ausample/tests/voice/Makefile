# Last edited on 2024-12-21 02:53:27 by stolfi

PROG := ausample
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# format:   FILEPREFIX:TAKE:EVERY:SPLICE

TESTS := \
  sony-recorder-sessao:000002.000:000180.000:00.200
  
OTHER_TESTS := \
  sony-recorder-stolfi:000001.500:000004.000:00.200

.PHONY:: loop run-single clean-single

ACTIONS := no-action

IN_DIR := in
OUT_DIR := out

all: 
	${MAKE} ACTIONS=single-run loop
        
clean::
	${MAKE} ACTIONS=single-clean loop

loop:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST="$${tt}" single-clean ; \
          ${MAKE} TEST="$${tt}" ${ACTIONS} ; \
        done

######################################################################
# Recursive make - caller must define ${TEST}

TEST := 
ifneq "/${TEST}" "/"

TESTARGS := ${subst :, ,${TEST}}
PREFIX   := ${word 1,${TESTARGS}}
TAKE     := ${word 2,${TESTARGS}}
EVERY    := ${word 3,${TESTARGS}}
SPLICE   := ${word 4,${TESTARGS}}

INP_AU_FILE  := ${IN_DIR}/${PREFIX}-inp.au
INP_OGG_FILE := ${IN_DIR}/${PREFIX}-inp.ogg
INP_PNG_FILE := ${IN_DIR}/${PREFIX}-inp.png

SMP_AU_FILE  := ${OUT_DIR}/${PREFIX}-smp.au
SMP_PNG_FILE := ${OUT_DIR}/${PREFIX}-smp.png
SMP_OGG_FILE := ${OUT_DIR}/${PREFIX}-smp.ogg

single-run: ${SMP_AU_FILE} ${SMP_PNG_FILE}

.SUFFIXES:

${SMP_AU_FILE}: ${INP_AU_FILE} ${PROGDIR}/${PROG} Makefile
	cat ${INP_AU_FILE} \
          | ${PROGDIR}/${PROG} \
              -take ${TAKE} -every ${EVERY} \
              -splice ${SPLICE} \
              -unit second \
          > ${SMP_AU_FILE}

#----------------------------------------------------------------------
# Sound file cleanup:

%-inp.au: %-inp.ogg Makefile
	sox -V $*-inp.ogg -e signed -b 16 $*-inp.au
	ls -l $*-inp.au

#----------------------------------------------------------------------
# Sound spectrograms:

%.png: %.au Makefile
	cat $*.au \
          | spegram \
              -input au \
              -window hann 1024 -stride 256 \
              -output pgm -range 1.0e-3 1.0e+1 \
          | pamflip -ccw \
          | pnminvert \
          | pnmnorm -bpercent 0 -wpercent 0 \
          | pnmtopng \
          > $*.png
	display $*.png

#----------------------------------------------------------------------
# Cleanup:

single-clean: 
	rm -f ${SMP_AU_FILE} ${SMP_OGG_FILE} ${SMP_PNG_FILE}

endif
# End ${TEST} section.
######################################################################
