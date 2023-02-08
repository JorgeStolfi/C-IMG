# Last edited on 2023-02-08 00:48:02 by stolfi

PROG := spegram
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

OGG_FILES := \
  data/snd-fausto-full.ogg \
  data/snd-sony-recorder-islene.ogg

AU_FILES :=  ${subst data/,out/, ${subst .ogg,.au,${OGG_FILES}}}

PGM_FILES := ${subst data/,out/, ${subst .ogg,.pgm,${OGG_FILES}}}

PNG_FILES := ${subst data/,out/, ${subst .ogg,.png,${OGG_FILES}}}

TXT_FILES := ${subst data/,out/, ${subst .ogg,.txt,${OGG_FILES}}}

# all: ${TXT_FILES}
all: ${PNG_FILES}

######################################################################
# OGG to AU conversion:

out/snd-%.au: data/snd-%.ogg 
	mkdir -p out
	sox -V data/snd-$*.ogg -e signed -b 16 out/snd-$*.au

######################################################################
# Sound spectrograms (image):

out/snd-%.png: ${PROGDIR}/${PROG} out/snd-%.au
	cat out/snd-$*.au \
          | ${PROGDIR}/${PROG} \
              -input au \
              -window hann 1024 -stride 256 \
              -output pgm -range 1.0e-3 1.0e+1 \
          | pnmflip -ccw \
          | pnminvert \
          | pnmnorm -bpercent 0 -wpercent 0 \
          | pnmtopng \
          > out/snd-$*.png
	display out/snd-$*.png

######################################################################
# Sound spectrograms (text):

out/snd-%.txt: ${PROGDIR}/${PROG} out/snd-%.au
	cat out/snd-$*.au \
          | ${PROGDIR}/${PROG} \
              -input au \
              -window hann 4096 -stride 1024 \
              -output ascii -range 0 1 \
          > out/snd-$*.txt
	plot_spectrogram.sh out/snd-$*.txt

######################################################################
# Cleanup:

clean::
	rm -f ${AU_FILES} ${PGM_FILES} ${PNG_FILES} ${TXT_FILES}
