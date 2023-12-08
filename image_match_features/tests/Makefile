# Last edited on 2020-11-06 20:11:29 by jstolfi

PROGDIR := ..
PROG = image_match_features

# Test runs
# Format: {IMG1}/{IMG2}/{HIX}/{HIY}/{DNAME}
# Where {IMG1} and {IMG2} are image names (without dir and ext), 
# {HIX,HIY} are the client coordinates of the upper corner of the image domain,
# and {DNAME} is the name of the input feature spec file(ditto).  

INDIR := in
OUTDIR := out

TESTS := \
  rgb-600-01/rgb-600-02/600/600/ft-600-01-02

ALL_TESTS := \
  rgb-128-01/rgb-128-02/1280/1280/ft-128-01-02 \
  gray-128-01/gray-128-02/1280/1280/ft-128-01-02 \
  rgb-128-01/rgb-128-03/1280/1280/ft-128-01-03

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

.PHONY:: all-tests single-test

all: check

check: clean all-tests single-test

all-tests:
	for tt in ${TESTS} ; do \
          ${MAKE} TEST="$$tt" single-test ; \
        done
        
########################################################################
# Section for recursive "make".
# Client must define ${TEST}

ifneq "/${TEST}" "/"

TEST_ARGS := ${subst /, ,${TEST}}

IMG1 :=   ${word 1, ${TEST_ARGS}}

IMG2 :=   ${word 2, ${TEST_ARGS}}

HIX :=    ${word 3, ${TEST_ARGS}}

HIY :=    ${word 4, ${TEST_ARGS}}

DNAME :=  ${word 5, ${TEST_ARGS}}

INFILE := ${INDIR}/${DNAME}.txt
OUTFILE := ${OUTDIR}/${DNAME}_adj.txt

FTRPREF := ${DNAME}
               
single-test: ${PROGDIR}/${PROG}
	rm -f ${OUTFILE} ${OUTDIR}/${FTRPREF}_*.png
	${PROGDIR}/${PROG} \
          -image1 ${INDIR}/${IMG1}.png 0 0 ${HIX} ${HIY} \
          -image2 ${INDIR}/${IMG2}.png 0 0 ${HIX} ${HIY} \
          -fix 1 \
          -maxIter 5 \
          -maxAdjust 0.75 \
          -writeFeatures ${FTRPREF} \
          -verbose \
          < ${INFILE} \
          > ${OUTFILE}
	cat ${OUTFILE}
	for img in ${OUTDIR}/${FTRPREF}_*.png ; do \
          display -filter box -resize '200%' $${img} ; \
        done

endif
# End ${TEST} section
########################################################################
