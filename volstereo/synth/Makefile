# Generates some synthetic stereo image pairs.
# Last edited on 2008-02-04 12:54:06 by stolfi

# Dimensions of mirr/021/r1-s.jpg
WIDTH := 512
HEIGHT := 384

# Same width as the "m12" image
# WIDTH := 349
# HEIGHT := 686

NRAYS = 2

OTHERINPUTS := \
  frag-stand.inc \
  frag-camera-whole.inc \
  frag-lights.inc \
  frag-monolith.inc \
  frag-frag.inc \
  frag-frag-grains.inc

# No need to change below here

NAME = main

POVPUB = ${STOLFIHOME}/programs/c-curitiba/povray-3.50c-1/PUB
POVRAY = ${POVPUB}/i686-Linux-2.2/bin/povray
POVINC = ${POVPUB}/share/povray-3.5/include
POVTTF = ${STOLFIHOME}/PUB/povray/tt-fonts
IMVIEW = /usr/X11R6/bin/display
CONVERT = /usr/X11R6/bin/convert

.PHONY:: clean  image showimage all-images index

all: image showimage

clean:: 
	-/bin/rm -f ${NAME}.png ${NAME}.gif

image: ${NAME}.gif

showimage: ${NAME}.gif
	-${IMVIEW} ${NAME}.gif

${NAME}.gif: ${NAME}.png
	-/bin/rm -f ${NAME}.gif
	${CONVERT} -gamma 2.2 ${NAME}.png ${NAME}.gif

${NAME}.png: ${NAME}.pov ${OTHERINPUTS}
	-/bin/rm -f ${NAME}.png
	nice ${POVRAY} \
            +FN +Q9 +MB1 \
            +W${WIDTH} +H${HEIGHT} \
            +AM1 +A0.0 +R${NRAYS} \
            +D +SP32 +EP4 \
            +L${POVINC} \
            +L${POVTTF} \
	    +I${NAME}.pov \
	    +O${NAME}.png \
	  3>&2 > ${NAME}.log

frag-frag-grains.inc: make-fragment-model
	make-fragment-model < /dev/null > frag-frag-grains.inc
        
