# Last edited on 2023-01-14 12:10:47 by stolfi

PROGDIR := ..
PROG := pnmspectrum

TESTS_A := abi2 ana2 chv1 eij1 lee1
TESTS_B := cpm2 hcg8 nie1

# TAG := A
# TESTS := ${TESTS_A}

TAG := B
TESTS := ${TESTS_B}

TABLE_TYPE := -exact

IN_IMAGES := ${addprefix data/, ${addsuffix .ppm,${TESTS}}}

OT_PREFIX := out/spectrum-${TAG}
OT_IMAGE := ${OT_PREFIX}.ppm
OT_TABLE := ${OT_PREFIX}.txt

all: ${PROGDIR}/${PROG} ${IN_IMAGES} out
	${PROGDIR}/${PROG} \
             -outputImage ${OT_IMAGE} \
               -center \
               -scale log 1.0e-6 \
             -outputTable ${OT_TABLE} ${TABLE_TYPE} \
             -verbose \
             ${IN_IMAGES}
	-display ${OT_IMAGE} ${IN_IMAGES}

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make
