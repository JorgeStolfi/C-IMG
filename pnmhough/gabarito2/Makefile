# Makefile pov-ray
# Last edited on 2008-02-04 12:52:15 by stolfi

WIDTH := 240
HEIGHT := 240

NRAYS = 5

OTHERINPUTS = 

# No need to change below here

NAME = main

ICONNAME = icon

POVPUB = ${STOLFIHOME}/programs/c-curitiba/povray-3.50c-1/PUB
POVRAY = ${POVPUB}/i686-Linux-2.2/bin/povray
POVINC = ${POVPUB}/share/povray-3.5/include
POVTTF = ${STOLFIHOME}/PUB/povray/tt-fonts
IMVIEW = /usr/X11R6/bin/display
CONVERT = /usr/X11R6/bin/convert

.PHONY:: image showimage all-images index

all: image showimage

clean::
	-/bin/rm -f ${NAME}.png ${NAME}.ppm

image: ${NAME}.ppm

showimage: ${NAME}.ppm
	-${IMVIEW} ${NAME}.ppm

${NAME}.ppm: ${NAME}-raw.ppm Makefile
	-/bin/rm -f ${NAME}.ppm
	${CONVERT} ${NAME}-raw.ppm -blur 5x0.5 -colorspace RGB ${NAME}.ppm

${NAME}-raw.ppm: ${NAME}.pov ${OTHERINPUTS}
	-/bin/rm -f ${NAME}.png
	nice ${POVRAY} \
            +FP +Q9 \
            +W${WIDTH} +H${HEIGHT} \
            +AM1 +A0.0 +R${NRAYS} \
            +D +SP32 +EP4 \
            +L${POVINC} \
            +L${POVTTF} \
	    +I${NAME}.pov \
	    +O${NAME}-raw.ppm \
	  3>&2 > ${NAME}.log

