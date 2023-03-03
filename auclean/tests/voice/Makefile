# Last edited on 2023-03-02 19:25:59 by stolfi

PROG := auclean
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

TESTS := \
  sony-recorder-islene
  
LATER_TESTS := \
  sony-recorder-rdahab \
  sony-recorder-stolfi \
  sony-recorder-gsouza \
  sony-recorder-islene \
  sony-recorder-sessao
  
.PHONY:: denoise make-pngs show-pngs
  
NOISE_FILE_1 := in/sony-recorder-noise1-raw.au
NOISE_FILE_2 := in/sony-recorder-noise2-raw.au
NOISE_FILE_3 := in/sony-recorder-noise3-raw.au
NOISE_FILE_4 := in/sony-recorder-noise4-raw.au

# Moment in file "islene" with voice:
# DEBUG_TIME := 12.66

# Moment in file "islene" with noise only:
DEBUG_TIME := 2.96

DEBUG_PREFIX := out/debug
DEBUG_FILE_01 := ${DEBUG_PREFIX}_01.pwr

RAW_AU_FILES := ${addprefix in/, ${addsuffix -raw.au,${TESTS}}}
RAW_PNG_FILES := ${addprefix in/, ${addsuffix -raw.png,${TESTS}}}

CLN_AU_FILES := ${addprefix out/, ${addsuffix -cln.au,${TESTS}}}
CLN_PNG_FILES := ${addprefix out/, ${addsuffix -cln.png,${TESTS}}}

NSE_AU_FILES := ${addprefix out/, ${addsuffix -nse.au,${TESTS}}}
NSE_PNG_FILES := ${addprefix out/, ${addsuffix -nse.png,${TESTS}}}

all: make-pngs show-pngs

.SUFFIXES:

######################################################################
# Sound file denoising

denoise: ${CLN_AU_FILES} ${NSE_AU_FILES}

out/%-cln.au out/%-nse.au: clean in/%-raw.au ${PROGDIR}/${PROG} Makefile plot_debug_pwr.sh
	cat in/$*-raw.au \
          | ${PROGDIR}/${PROG} \
              -window 2048 -overlap 2 \
              -noiseFile ${NOISE_FILE_1} \
              -noiseFile ${NOISE_FILE_2} \
              -noiseFile ${NOISE_FILE_3} \
              -noiseFile ${NOISE_FILE_4} \
              -debugTime ${DEBUG_TIME} \
              -writeChaff \
              -outPrefix out/$*
	if [[ -s ${DEBUG_FILE_01} ]]; then \
          for dbgf in ${DEBUG_PREFIX}_??.pwr ; do \
            ./plot_debug_pwr.sh $$dbgf 8000 ; \
          done ; \
        fi

NO_OPTIONS := \
              -kill 0.0 254 100 \
              -fmax 3900 \
              -fmin   50 \
              -noiseInput          

######################################################################
# Sound spectrograms:

make-pngs: ${CLN_PNG_FILES} ${NSE_PNG_FILES}

show-pngs: make-pngs
	for fot in ${CLN_PNG_FILES} ; do \
          fin="$${fot/out/in}"; fin="$${fin/cln/raw}"; \
          fns="$${fot/cln/nse}"; \
          display -loop 0 -delay 100x1 $${fin} $${fot} $${fns} ; \
        done

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

######################################################################
# Cleanup:

clean::
	rm -f \
          ${CLN_AU_FILES} ${CLN_PNG_FILES} \
          ${NSE_AU_FILES} ${NSE_PNG_FILES} \
          ${DEBUG_PREFIX}_??.pwr
