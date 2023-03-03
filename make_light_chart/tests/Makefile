# Last edited on 2023-03-02 22:35:55 by stolfi

PROG := make_light_chart
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make


PSVIEW := atril

OUT_PREFIX := out/target

all:
	${PROG} -radius 5.0 -rings 4 -outPrefix ${OUT_PREFIX}
	${PSVIEW} ${OUT_PREFIX}.eps        
	${PROG} -radius 10.0 -rings 4 -outPrefix ${OUT_PREFIX}-n -showSpotNumbers
	${PSVIEW} ${OUT_PREFIX}-n.eps
