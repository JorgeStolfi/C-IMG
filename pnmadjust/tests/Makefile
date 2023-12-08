# Last edited on 2010-07-02 17:06:40 by stolfi

PROG := pnmadjust
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  1999-stmalo-02

TESTS_LATER :=
  
TESTSBUG := \
  ext-noronha-1 \
  oliva-1
  
NONTESTS := \
  gray-50

UNIFTESTS_OK := \
  1999-stmalo-02 \
  abistolfi-2 \
  ariadne-13 \
  arnaldo-4 \
  beatriz-3 \
  candida-1 \
  cmrubira-1 \
  eliane-1 \
  ext-avg-1 \
  hcgl-12 \
  hcgl-16 \
  keidi-6 \
  lhf-1 \
  ragazzi-4 \
  stolfi-4 \
  wolsey-1 \
  iclabs-08
  
VARTESTS_OK := \
  lucchesi-3 \
  2000-talossa-03 \
  cecilia-4 \
  ext-jmason-1 \
  freitas-1 \
  gray-50 \
  marta-3 \
  pannain-1

VARTESTS_BUG :=

.PHONY:: clean single single-clean single-img

all:
	for t in ${TESTS} ; do \
          ${MAKE} TEST=$$t single ; \
        done

clean::
	rm -f *-o.* *-c.*

######################################################################
# Recursive make - must define ${TEST}

TEST := TEST.IS.UNDEFINED
ifneq "/${TEST}" "/TEST.IS.UNDEFINED"

PPMRAW := data/${TEST}.ppm.gz
JPGOLD := data/${TEST}-o.jpg

PPMCOR := out/${TEST}-c.ppm.gz
JPGCOR := out/${TEST}-c.jpg

single: single-clean single-img

single-clean:
	rm -f ${PPMCOR} ${JPGCOR}

single-img: ${PPMCOR} ${JPGCOR} ${JPGOLD}
	display -title '%f' ${PPMCOR} ${PPMRAW} ${JPGCOR} ${JPGOLD} 

${PPMCOR}: ${PROGDIR}/${PROG}
	cat ${PPMRAW} \
          | gunzip \
          | ${PROGDIR}/${PROG} \
              `cat ${TEST}.parms | sed -e 's/[#].*//g'` \
          | gzip \
          > ${PPMCOR}

%-c.jpg: %-c.ppm.gz
	cat $< \
          | gunzip \
          | pnmgamma 2.2 \
          | convert -quality 98 PPM:- $*-c.jpg

%-o.jpg: %.ppm.gz
	cat $< \
          | gunzip \
          | pnmgamma 2.2 \
          | convert -quality 98 PPM:- $*-o.jpg
        
endif
# end ${TEST}
######################################################################

######################################################################
# Rebuilding the program:

${PROGDIR}/${PROG}: ${PROGDIR} ${PROGDIR}/${PROG}.c
	cd ${PROGDIR} && ${MAKE} ${PROG}
