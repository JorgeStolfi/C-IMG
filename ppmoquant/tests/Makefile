# Last edited on 2009-01-07 02:19:32 by stolfi

PROG := ppmoquant
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

CMAP := ${STOLFIHOME}/PUB/colormaps/standard-6x6x6-colormap-gm10.ppm

all: ${PROGDIR}/${PROG} 
	${PROGDIR}/${PROG} \
            -fs -map ${CMAP} \
            -transparent 64 64 64 \
            black.ppm white.ppm \
            > trans.ppm
	cat trans.ppm \
          | pnmgamma 1.8 \
          | ppmtogif -transparent rgb:76/76/76 \
          > trans.gif
clean::

