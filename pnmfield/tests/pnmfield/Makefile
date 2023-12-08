# Last edited on 2009-01-06 17:51:20 by stolfi

PROG := pnmfield
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  air-canyon-1-black.ppm \
  air-canyon-1-white.ppm

TESTS_NOT := \
  rgb-lin-ramp.ppm \
  rgb-lin-uniform.ppm \
  rgb-lin-hWave.ppm \
  rgb-lin-vWave.ppm \
  rgb-lin-wave.ppm \
  rgb-lin-wavePair.ppm \
  gry-lin-wave.ppm \
  gry-lin-ramp.ppm \
  rgb-log-wavePair.ppm \
  rgb-log-ramp.ppm \
  rgb-lin-wave-igam.ppm \
  rgb-lin-wave-ogam.ppm

all: ${TESTS}

clean::
	rm -f ${TESTS}

SHOW := display
HRULE := echo ' --------------------------------------------- '
        
rgb-lin-ramp.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field ramp \
	      20 30  127 127 255 / 255 \
	      50 60  127 255 127 / 255 \
	      70 20  255 127 127 / 255 \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-lin-uniform.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field uniform \
	      255 127 000 / 255 \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-lin-hWave.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field hWave \
	      10  191 032 063 / 255 \
	      20  063 223 191 / 255 \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-lin-vWave.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field vWave \
	      10  191 032 063 / 255 \
	      20  063 223 191 / 255 \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-lin-wave.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field wave \
	      10 20  191 032 063 / 255 \
	      20 30  063 223 191 / 255 \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-lin-wavePair.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field wavePair \
	      10 20  255 255 255 / 255 \
	      20 40  255 191 127 / 255 \
	      30 10  191 255 127 / 255 \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

gry-lin-wave.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field wave \
	      10 20  191 032 063 / 255 \
	      20 30  063 223 191 / 255 \
	    -gray \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

gry-lin-ramp.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field ramp \
	      20 30  000 000 000 / 255 \
	      30 40  127 000 127 / 255 \
	      60 20  000 127 127 / 255 \
	    -gray \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-log-wavePair.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field wavePair \
	      10 20  255 255 255 / 255 \
	      20 40  255 191 127 / 255 \
	      30 10  191 255 127 / 255 \
	    -log \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-log-ramp.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field ramp \
	      20 30  127 127 255 / 255 \
	      50 60  127 255 127 / 255 \
	      70 20  255 127 127 / 255 \
	    -log \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-lin-wave-igam.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
            -inGamma 2.0 2.0 2.0 \
            -field wave \
	      10 20  255 255 255 / 255 \
	      20 30  063 127 191 / 255 \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

rgb-lin-wave-ogam.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field wave \
	      10 20  063 127 191 / 255 \
	      20 30  255 255 255 / 255 \
            -outGamma 2.0 2.0 2.0 \
	    80 80 \
	  > $@
	${SHOW} $@
	@${HRULE}

air-canyon-1-black.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field wavePair \
              0030 0220   039 062 094 / 255 \
              0030 0000   108 152 176 / 255 \
              0260 0220   070 109 133 / 255 \
           490 289 \
	  > $@
	${SHOW} $@
	@${HRULE}

air-canyon-1-white.ppm: ${PROGDIR}/${PROG} Makefile
	${PROGDIR}/${PROG} \
	    -field wavePair \
              0220 0020   243 227 226 / 255 \
              0000 0220   201 207 206 / 255 \
              0200 0260   221 205 188 / 255 \
	    490 289 \
	  > $@
	${SHOW} $@
	@${HRULE}

# 
# 	    
#   0034 0228  045 071 102 / 255
#   0312 0254  074 114 134 / 255
#   0448 0224  059 091 119 / 255
#   0020 0098  094 123 142 / 255
#   0180 0101  121 159 174 / 255
#   0432 0138  091 121 145 / 255
#   0017 0022  107 142 166 / 255
#   0269 0072  144 177 184 / 255
#   0436 0012  170 195 199 / 255
