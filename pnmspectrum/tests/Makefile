# Last edited on 2024-10-17 12:28:01 by stolfi

PROGDIR := ..
PROG := pnmspectrum

.PHONY:: all single-test

TESTS_A := photos-people:abi2.ppm+ana2.ppm+chv1.ppm+eij1.ppm+lee1.ppm
TESTS_B := photos-stripescpm2.ppm+hcg8.ppm+nie1.ppm
TESTS_C := 
TESTS_T := texto-double.ppm+texto-shrink.ppm+texto-stretch.ppm+texto-normal.ppm

# TESTS := photos-people:abi2.ppm+ana2.ppm+chv1.ppm+eij1.ppm+lee1.ppm
# TESTS := photos-stripescpm2.ppm+hcg8.ppm+nie1.ppm

TESTS := grittie:fimt-ns3/grittie.pgm chopsea:fimt-ns3/chopsea.pgm bullsqr:fimt-ns3/bullsqr.pgm bullsex:fimt-ns3/bullsex.pgm

# TESTS := grittie-sma:fimt-ns3/grittie-sma.ppm chopsea-sma:fimt-ns3/chopsea-sma.ppm bullsqr-sma:fimt-ns3/bullsqr-sma.ppm bullsex-sma:fimt-ns3/bullsex-sma.ppm
  
# TESTS := texts:texto-double.ppm+texto-shrink.ppm+texto-stretch.ppm+texto-normal.ppm

TABLE_TYPE := -exact

all:
	for test in ${TESTS}; do \
          ${MAKE} TEST="$${test}" single-test ; \
        done

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

# ######################################################################
# Recursive make section for ${TEST}

ifneq "/${TEST}" "/"

TEST_ARGS := ${subst :, ,${TEST}}
TAG :=     ${word 1, ${TEST_ARGS}}
IMGLIST := ${word 2, ${TEST_ARGS}}

IN_IMAGES := ${addprefix data/, ${subst +, ,${IMGLIST}}}

OT_PREFIX := out/spectrum-${TAG}
OT_IMAGE := ${OT_PREFIX}.ppm
OT_TABLE := ${OT_PREFIX}.txt

single-test: ${PROGDIR}/${PROG} ${IN_IMAGES} out
	${PROGDIR}/${PROG} \
             -outputImage ${OT_IMAGE} \
               -center \
               -scale log 1.0e-6 \
             -outputTable ${OT_TABLE} ${TABLE_TYPE} \
             -verbose \
             ${IN_IMAGES}
	-display -title '%f' -filter Box -resize '100%' ${OT_IMAGE} ${IN_IMAGES}

endif

# end ${TEST} section
# ######################################################################
