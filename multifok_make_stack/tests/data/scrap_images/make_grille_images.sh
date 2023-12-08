#! /bin/bash
# Last edited on 2018-07-02 08:25:35 by stolfilocal

# Create a field of small crosses and dots, using the floor grid
# created for the microscope UFF project.

#---------------------------------------------------------------------------
# Getting the source images
#---------------------------------------------------------------------------

# Get the PNG files of the test grid sheets:
gdir="programs/c/IMG/make_grid_sheet/tests/out"

ifiles=( \
    grid_54x54_2052x2052.png \
    grid_64x64_2048x2048.png \
    grid_58x58_2088x2088.png \
  )
 
for ifile in ${ifiles[@]} ; do 
    cp -v ~/${gdir}/${ifile} ${ifile}
    identify ${ifile}
    display ${ifile}
  done

# Checking for symmetry and accurate positioning of marks:
for ff in ${ifiles[@]}; do
  convert ${ff} .foo.pgm
  pamflip -rotate180 .foo.pgm > .bar.pgm
  pnmxarith -difference .foo.pgm .bar.pgm | pgmnorm -bsingle -wsingle > .diff.pgm
  display .diff.pgm
done
# The images are totally symmetric.  

#---------------------------------------------------------------------------
# Creating derived images
#---------------------------------------------------------------------------

# Cropping a an image with 27 x 27 cells at 38 pixels per cell. 
ifile="grid_54x54_2052x2052.png"
ofile1="grid_27x27_1026x1026.png"
convert ${ifile} -crop '1026x1026+513+513' ${ofile1}
identify ${ofile1}
display ${ofile1}

# Create an image with 27 x 27 of 19 pixels per cell. 
ifile="grid_27x27_1026x1026.png"
ofile2="grid_27x27_513x513.png"
convert ${ifile} -resize '513x513' ${ofile2}
identify ${ofile2}
display ${ofile2}

# Create an image with 57 x 57 cells of 9 pixels per cell. 
ifile="grid_58x58_2088x2088.png"
ofile3="grid_57x57_513x513.png"
convert ${ifile} -crop '2052x2052+18+18' -resize '513x513' ${ofile3}
identify ${ofile3}
display ${ofile3}

