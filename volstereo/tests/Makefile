# Last edited on 2017-06-23 23:47:46 by stolfilocal

PROGS := \
  draw_corresp_pairs.gawk \
  make_coincidence_map.gawk

PROGDIR := ..

# include ${STOLFIHOME}/programs/c/GENERIC-SCRIPTS-TEST.make

VARWIDTH1 := \
  244:1 \
  244:3 \
  244:5 \
  244:7 \
  244:9
  
VARWIDTH2 := \
  140:1 \
  140:3 \
  140:5 \
  140:7 \
  140:9
  
MISCFRAMES1 := \
  244:1 \
  244:9

MISCFRAMES2 := \
  090:3 \
  095:3 \
  153:3 \
  244:3 \
  244:5 \
  244:7 \
  140:1

VARCOL1 := \
  236:7 \
  238:7 \
  240:7 \
  242:7 \
  246:7 \
  248:7 \
  250:7 \
  252:7 \
  254:7

VARCOL2 := \
  040:7 \
  045:7 \
  050:7 \
  055:7 \
  060:7 \
  065:7 \
  070:7 \
  075:7 \
  080:7 \
  085:7 \
  090:7 \
  095:7 \
  100:7 \
  105:7 \
  110:7 \
  115:7 \
  120:7 \
  125:7 \
  130:7 \
  135:7 \
  140:7 \
  145:7 \
  150:7

# FRAMES := ${VARCOL1}
# FRAMES := ${VARCOL2}
FRAMES := ${MISCFRAMES1}
# FRAMES := ${VARWIDTH2}

MROW := 515
ROWS := 200

SHOWFRAMES := 1
MAKEMOVIE := 0

all: frames movie

frames:
	rm -f .frames
	for t in ${FRAMES}; do \
          ${MAKE} FRAMESPECS=$$t single ; \
        done

clean-all:
	rm -f out/*.pcmds out/*-cmap.* out/*-movie*.gif 
	rm -f out/*-???-???-?.ppm

######################################################################
# Recursive make - must define ${FRAMESPECS} = {MCOL}:{COLS}
# = middle column and number of columns.

FRAMESPECS := FRAMESPECS.IS.UNDEFINED
ifneq "/${FRAMESPECS}" "/FRAMESPECS.IS.UNDEFINED"

SPECFIELDS := ${subst :, ,${FRAMESPECS}}
MCOL := ${word 1,${SPECFIELDS}}
COLS := ${word 2,${SPECFIELDS}}

# Initial column and row:
FCOL := ${shell echo ${MCOL} | gawk -v m=${MCOL} -v n=${COLS} 'BEGIN{print m-int(n/2);}'}
FROW := ${shell echo ${MROW} | gawk -v m=${MROW} -v n=${ROWS} 'BEGIN{print m-int(n/2);}'}

TAG := ${MCOL}-${MROW}-${COLS}

IMG_PPM         := in/m21.ppm
IMG_PAIRS       := in/m21.pairs
IMG_SLICED_PPM  := out/m21-${TAG}.ppm

CMAP_ANNOT_CMDS := out/m21-${TAG}.pcmds
CMAP_PGM        := out/m21-${TAG}-cmap.pgm
CMAP_ANNOTATED  := out/m21-${TAG}-cmap.png

MOVIE_GIF       := out/m21-${MROW}-movie.gif

SCALE := 2

# NORM := pgmlocnorm -radius 7
NORM := pgmnorm -bvalue 0 -wvalue 2000
# NORM := cat

REDUCE := ppmtopgm | pnmnoraw
# REDUCE := cat

# single: ; echo "${MCOL} ${COLS} ${FCOL}  ${MROW} ${ROWS} ${FROW}"
single: ${CMAP_ANNOTATED}

${CMAP_ANNOTATED}: ${CMAP_PGM} ${CMAP_ANNOT_CMDS} Makefile
	cat ${CMAP_PGM} \
          | pnminvert \
          | convert PGM:- \
              -fill yellow \
              -draw '@'${CMAP_ANNOT_CMDS} \
              PNG:- \
          > ${CMAP_ANNOTATED}
	if [ ${SHOWFRAMES} ] ; then\
          display ${CMAP_ANNOTATED} ${CMAP_PGM}; \
        fi

${CMAP_PGM}: ${IMG_SLICED_PPM} ${PROGDIR}/make_coincidence_map.gawk Makefile
	cat ${IMG_SLICED_PPM} \
          | ${REDUCE} \
          | ${PROGDIR}/make_coincidence_map.gawk \
          | ${NORM} \
          > ${CMAP_PGM}
	echo ${CMAP_PGM} >> .frames

${IMG_SLICED_PPM}: ${IMG_PPM} Makefile
	cat ${IMG_PPM} \
          | pnmcut ${FCOL} ${FROW} ${COLS} ${ROWS} \
          | pnmscale -yscale ${SCALE} \
          | pnmnoraw \
          > ${IMG_SLICED_PPM}

${CMAP_ANNOT_CMDS}: ${IMG_PAIRS} ${PROGDIR}/draw_corresp_pairs.gawk Makefile
	cat ${IMG_PAIRS} \
          | ${PROGDIR}/draw_corresp_pairs.gawk \
              -v frow=${FROW} -v rows=${ROWS} \
              -v fcol=${FCOL} -v cols=${COLS} \
              -v scale=${SCALE} \
          > ${CMAP_ANNOT_CMDS}

endif
# end ${FRAMESPECS}
######################################################################

movie: .frames
	rm -f ${MOVIE_GIF}
	if [ ${MAKEMOVIE} ]; then \
          convert \
            -delay 40 -loop 10 \
            `cat .frames` \
            ${MOVIE_GIF} ; \
        fi
	if [ ${SHOWMOVIE} ] ; then\
          display -delay 50 ${MOVIE_GIF}; \
        fi
