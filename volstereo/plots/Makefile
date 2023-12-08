# generates coincidence plots from a stereo image pair
# Last edited on 2023-02-04 07:06:42 by stolfi

OTHERPLOTS := \
  140:433:502:533:602

NONPLOTS := \
  244:465:525:535:605

PLOTS := ${OTHERPLOTS}

all: plots

plots:
	rm -f .plots
	for t in ${PLOTS}; do \
          ${MAKE} PLOTSPECS=$$t single ; \
        done

######################################################################
# Recursive make - must define
# ${PLOTSPECS} = {COL}:{IROW0}:{FROW0}:{IROW1}:{FROW1}
# where ${IROW0}..${FROW0} is paired to {IROW1}..{FROW1} reversed.

PLOTSPECS := PLOTSPECS.IS.UNDEFINED
ifneq "/${PLOTSPECS}" "/PLOTSPECS.IS.UNDEFINED"

SPECFIELDS := ${subst :, ,${PLOTSPECS}}
COL := ${word 1,${SPECFIELDS}}
IROW0 := ${word 2,${SPECFIELDS}}
FROW0 := ${word 3,${SPECFIELDS}}
IROW1 := ${word 4,${SPECFIELDS}}
FROW1 := ${word 5,${SPECFIELDS}}

ROWS0 := ${shell echo "$$((${FROW0} - ${IROW0}))"}
ROWS1 := ${shell echo "$$((${FROW1} - ${IROW1}))"}

IMG_NAME        := m21

IMG_PPM         := ${IMG_NAME}.ppm
IMG_PAIRS       := ${IMG_NAME}.pairs

PLOT_EPS        := ${IMG_NAME}-${COL}.eps
PLOT_CMDS       := ${IMG_NAME}-${COL}.plot
PLOT_DATA_0     := ${IMG_NAME}-${COL}-0.data
PLOT_DATA_1     := ${IMG_NAME}-${COL}-1.data

REDUCE := ppmtopgm
CHANNELS = 1

# REDUCE := cat
# CHANNELS := 3

# single: ; echo "${COL} ${IROW0}..${FROW0}(${ROWS0}) ${IROW1}..${FROW1}(${ROWS0})"
single: ${PLOT_EPS}

PSVIEW := evince

${PLOT_EPS}: ${PLOT_DATA_0} ${PLOT_DATA_1} ${PLOT_CMDS} Makefile
	gnuplot < ${PLOT_CMDS}
	${PSVIEW} ${PLOT_EPS}

${PLOT_CMDS}: make-plot-cmds Makefile
	make-plot-cmds \
            -v img="${IMG_NAME}"  -v col="${COL}" \
            -v irow0="${IROW0}" -v frow0="${FROW0}" \
            -v irow1="${IROW1}" -v frow1="${FROW1}" \
            -v channels=${CHANNELS} \
          < /dev/null \
          > ${PLOT_CMDS}
          
${PLOT_DATA_0}: make-plot-data ${IMG_PPM} Makefile
	pnmcut ${COL} ${IROW0} 1 ${ROWS0} ${IMG_PPM} \
          | ${REDUCE} | pnmnoraw \
          | make-plot-data \
          > ${PLOT_DATA_0}

${PLOT_DATA_1}: make-plot-data ${IMG_PPM} Makefile
	pnmcut ${COL} ${IROW1} 1 ${ROWS1} ${IMG_PPM} \
          | ${REDUCE} | pnmnoraw \
          | make-plot-data \
          > ${PLOT_DATA_1}

endif

movie: .plots
	rm -f ${IMG_NAME}-${MROW}-movie.gif
	if [ ${MAKEMOVIE} ]; then \
          convert \
            -delay 40 -loop 10 \
            `cat .plots` \
            ${IMG_NAME}-${MROW}-movie.gif ; \
        fi
