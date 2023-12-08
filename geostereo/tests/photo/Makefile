# Last edited on 2008-01-14 20:23:13 by stolfi

PROG := geostereo
PROGDIR := ../..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

all: fa-disp.pgm

clean::

fa-disp.pgm: \
          ${PROGDIR}/${PROG} \
          Makefile \
          f0a-e.pgm f1a-e.pgm f0a-s.pgm f1a-s.pgm
	${PROGDIR}/${PROG} \
	    -window 3 3 -nscales 0 \
	    -mindisp 0 -maxdisp +10 \
	    -output fa \
	    f0a-e.pgm f1a-e.pgm
	pnmxarith -subtract -offset 128 f0a-e.pgm f1a-e.pgm > fa-diff.pgm
	display \
          fa-disp.pgm fa-score.pgm \
          f0a-e.pgm f1a-e.pgm \
          fa-diff.pgm \
          f0a-s.pgm f1a-s.pgm 

%-e.pgm: %.jpg %.epmap Makefile
	convert $*.jpg PPM:- \
          | ppmtopgm \
          | pnmgtran \
	      -matrix `head -3 $*.epmap` \
	  | pnmscale 0.5 \
          | pgmnorm -bpercent 0 -wpercent 0 \
          > $*-e.pgm
	display $*-e.pgm $*.jpg
         
%-s.pgm: %.jpg
	convert $*.jpg PPM:- \
          | ppmtopgm \
          | pnmscale 0.25 \
          | pgmnorm -bpercent 0 -wpercent 0 \
          > $*-s.pgm
	display $*-s.pgm $*.jpg
         
%.epmap: %.grmap  compute-epi-map r3x3.gawk
	compute-epi-map -f r3x3.gawk \
	    -v width=2048 \
	    -v check=1 \
	  < $*.grmap > $*.epmap
	 
%.grmap: %.grpts  compute-grid-map r3x3.gawk
	compute-grid-map -f r3x3.gawk \
	    -v check=1 \
	  < $*.grpts > $*.grmap
