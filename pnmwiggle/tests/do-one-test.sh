#! /bin/bash
# Last edited on 2013-02-12 04:49:44 by stolfilocal

file="$1"; shift;
name=${file%%.*}
ext=${file##*.}
# echo "name = [${name}] ext = [${ext}]"
../pnmwiggle 0.25 in/${file} > out/${name}-025.${ext}
../pnmwiggle 1.00 in/${file} > out/${name}-100.${ext}
display out/${name}-*.${ext} in/${file}
