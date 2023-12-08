#! /bin/bash
# Last edited on 2018-07-02 09:28:05 by stolfilocal

# Creates a synthetic grid-like image and a height map
# with mounds of various shapes, smooth and discontinuous.

cp -v ../scrap_images/grid_57x57_513x513.png  main.png
convert ../scrap_images/ht_sevens_01.pgm      height.png
