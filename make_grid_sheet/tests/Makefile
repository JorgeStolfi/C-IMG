# Last edited on 2020-11-29 21:29:07 by jstolfi

PROG := make_grid_sheet
PROGDIR := ..

include ${STOLFIHOME}/programs/c/GENERIC-PROGS-TEST.make

PDFVIEW := evince
PNGVIEW := display

OUTFILES := \
  out/grid_54x54.eps \
  out/grid_54x54_2052x2052.png \
  out/grid_64x64_2048x2048.png \
  out/grid_58x58_2088x2088.png 

all: ${OUTFILES}

out/grid_54x54.eps: ${PROGDIR}/${PROG} 
	${PROGDIR}/${PROG} \
            -gridSize 54 54 \
            -markSpacing 1.0 \
            -markRadius 0.400 \
            -lineWidth 0.15 \
            -printerDPI 600 \
            -outPrefix out/grid_54x54
	${PDFVIEW} out/grid_54x54.eps        

out/grid_54x54_2052x2052.png: ${PROGDIR}/${PROG} 
	${PROGDIR}/${PROG} \
            -gridSize 54 54 \
            -markSpacing 1.0 \
            -markRadius 0.400 \
            -lineWidth 0.20 \
            -imagePXPMM 38 \
            -outPrefix out/grid_54x54_2052x2052
	identify out/grid_54x54_2052x2052.png
	display out/grid_54x54_2052x2052.png

out/grid_64x64_2048x2048.png: ${PROGDIR}/${PROG} 
	${PROGDIR}/${PROG} \
            -gridSize 64 64 \
            -markSpacing 1.0 \
            -markRadius 0.400 \
            -lineWidth 0.20 \
            -imagePXPMM 32 \
            -outPrefix out/grid_64x64_2048x2048
	identify out/grid_64x64_2048x2048.png
	display out/grid_64x64_2048x2048.png

out/grid_58x58_2088x2088.png: ${PROGDIR}/${PROG} 
	${PROGDIR}/${PROG} \
            -gridSize 58 58 \
            -markSpacing 1.0 \
            -markRadius 0.400 \
            -lineWidth 0.15 \
            -imagePXPMM 36 \
            -outPrefix out/grid_58x58_2088x2088
	identify out/grid_58x58_2088x2088.png
	display out/grid_58x58_2088x2088.png

# 0.16933333333333333333
# 0.08466666666666666666
# 0.04233333333333333333
# 0.08466666666666666666
