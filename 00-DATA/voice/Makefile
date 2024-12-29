# Last edited on 2023-03-02 19:24:55 by stolfi

RAWFILES := \
  sony-recorder-gsouza-raw.au \
  sony-recorder-islene-raw.au \
  sony-recorder-rdahab-raw.au \
  sony-recorder-sessao-raw.au \
  sony-recorder-stolfi-raw.au

NOISEFILES := \
  sony-recorder-noise1-raw.au \
  sony-recorder-noise2-raw.au \
  sony-recorder-noise3-raw.au \
  sony-recorder-noise4-raw.au

AUFILES := ${RAWFILES} ${NOISEFILES}

PNGFILES := ${subst .au,.png,${AUFILES}}

all: ${AUFILES} ${PNGFILES}

sony-recorder-noise1-raw.au: sony-recorder-islene-raw.au 
	extract_au_clip 420 520 2427 \
          < sony-recorder-islene-raw.au \
          > sony-recorder-noise1-raw.au

sony-recorder-noise2-raw.au: sony-recorder-islene-raw.au 
	extract_au_clip 940 1040 2427 \
          < sony-recorder-islene-raw.au \
          > sony-recorder-noise2-raw.au

sony-recorder-noise3-raw.au: sony-recorder-islene-raw.au 
	extract_au_clip 1274 1327 2427 \
          < sony-recorder-islene-raw.au \
          > sony-recorder-noise3-raw.au

sony-recorder-noise4-raw.au: sony-recorder-islene-raw.au 
	extract_au_clip 1645 1665 2427 \
          < sony-recorder-islene-raw.au \
          > sony-recorder-noise4-raw.au

%.png: %.au Makefile
	cat $*.au \
          | spegram \
              -input au \
              -window hann 1024 -stride 256 \
              -output pgm -range 1.0e-3 1.0e+1 \
          | pnmflip -ccw \
          | pnminvert \
          | pnmnorm -bpercent 0 -wpercent 0 \
          | pnmtopng \
          > $*.png
	display $*.png

clean: 
	rm -fv ${PNGFILES} ${NOISEFILES}
