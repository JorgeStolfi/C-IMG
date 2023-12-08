# Last edited on 2011-06-05 23:57:06 by stolfi

PROG := ppmtoymn
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TEST := app

ALL_TESTS := \
  grd \
  amg \
  agg \
  agh \
  app
  
# grd = grade de cores 6x6x6, com valores k*55 em cada calal, k=0..5
# amg = retrato, acredita-se gama = 1.0
# agg = retrato, acredita-se gama = 2.0
# agg = retrato, acredita-se gama = 0.5
# app = maçãs, gama desconhecido

IN_PPM := in-${TEST}.ppm

TARGETS := \
  out/${TEST}-ident.ppm \
  out/${TEST}-ymn.ppm \
  out/${TEST}-y.ppm \
  out/${TEST}-uv.ppm \
  out/${TEST}-u.ppm \
  out/${TEST}-v.ppm \
  out/${TEST}-clr.ppm \
  out/${TEST}-ny.ppm \
  out/${TEST}-nuv.ppm
  
MORE_TARGETS := \
  out/${TEST}-bri.ppm

.PHONY:: targets

all: clean targets
  

clean::
	rm -f ${TARGETS}
        
targets: ${TARGETS}
	display -title '%f' ${IN_PPM} ${TARGETS}


out/${TEST}-ident.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale +1.0 +1.0 +1.0 \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -inverse \
          > $@

out/${TEST}-ymn.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale +1.0 +1.0 +1.0 \
            ${IN_PPM} \
          > $@

out/${TEST}-y.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale +1.0 00.0 00.0 \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -inverse \
          > $@

out/${TEST}-uv.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale 00.0 +1.0 +1.0 \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -inverse \
          > $@

out/${TEST}-u.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale 00.0 +1.0 00.0 \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -inverse \
          > $@

out/${TEST}-v.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale 00.0 00.0 +1.0 \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -inverse \
          > $@

out/${TEST}-clr.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale +1.0 +1.5 +1.5 \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -inverse \
          > $@
          
out/${TEST}-ny.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale -1.0 +1.0 +1.0 \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -inverse \
          > $@

out/${TEST}-nuv.ppm: ${IN_PPM} ${PROGDIR}/${PROG}
	cat ${IN_PPM} \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -rgbscale +1.0 +1.0 +1.0 \
            -yuvscale +1.0 -1.0 -1.0 \
          | ${PROGDIR}/${PROG} \
            ${BIAS} \
            -inverse \
          > $@

          
