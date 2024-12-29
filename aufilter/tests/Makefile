# Last edited on 2024-12-20 19:42:09 by stolfi

PROG := aufilter
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  sony-recorder-islene
  
LATER_TESTS :=
  
NOISE_FILE := in/sony-recorder-noises-raw.au

# DEBUG_TIME := 12.66
DEBUG_TIME := 2.96
DEBUG_FILE := out/debug.pwr

RAW_AU_FILES := ${addprefix in/, ${addsuffix -raw.au,${TESTS}}}
RAW_PNG_FILES := ${addprefix in/, ${addsuffix -raw.png,${TESTS}}}

CLN_AU_FILES := ${addprefix out/, ${addsuffix -cln.au,${TESTS}}}
CLN_PNG_FILES := ${addprefix out/, ${addsuffix -cln.png,${TESTS}}}

# all: ${RAW_AU_FILES} ${RAW_PNG_FILES} ${CLN_AU_FILES} ${CLN_PNG_FILES}
all: ${RAW_AU_FILES} ${CLN_AU_FILES} ${CLN_PNG_FILES}

.SUFFIXES:

######################################################################
# Sound file denoiseup:

out/%-cln.au: in/%-raw.au ${PROGDIR}/${PROG} Makefile plot_debug_pwr.sh
	rm -f out/$*-cln.au out/$*-cln.png ${DEBUG_FILE}
	cat in/$*-raw.au \
          | ${PROGDIR}/${PROG} \
              -window 2048 -overlap 2 \
              -kill 0.0 254 100 \
              -fmax 3900 \
              -fmin   50 \
              -debugTime ${DEBUG_TIME} \
          > out/$*-cln.au
	if [[ -s ${DEBUG_FILE} ]]; then \
          ./plot_debug_pwr.sh ${DEBUG_FILE} ; \
        fi

NO_OPTIONS := \
              -noiseInput          

######################################################################
# Sound spectrograms:

out/%.png: out/%.au Makefile
	cat out/$*.au \
          | spegram \
              -input au \
              -window hann 1024 -stride 256 \
              -output pgm -range 1.0e-3 1.0e+1 \
          | pnmflip -ccw \
          | pnminvert \
          | pnmnorm -bpercent 0 -wpercent 0 \
          | pnmtopng \
          > out/$*.png
	fot="$*"; eom in/$${fot/cln/raw}.png out/$*.png

######################################################################
# Cleanup:

clean::
	rm -f ${CLN_AU_FILES} ${CLN_PNG_FILES} ${DEBUG_FILE}
