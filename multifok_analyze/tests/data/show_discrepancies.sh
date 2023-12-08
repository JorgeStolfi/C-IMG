#! /bin/bash 
# Last edited on 2018-01-04 00:24:05 by stolfilocal

# USAGE: ( cd {STACKDIR} ; show_discrepancies.sh )

function prepare_file() {
  ifile="$1"; shift;
  pfile="$1"; shift;
  echo "${ifile} --> ${pfile}"
  convert ${ifile} PGM:- \
    | pnmwfilter \
        -filter normalize 3.0 \
        -weights weights_5x5.pgm \
        -noise 0.08 \
        -maxval 65535 \
        -replicate \
    > ${pfile}
}

function compare_files() {
  fpgm="$1"; shift;
  rpgm="$1"; shift;
  opgm="$1"; shift;
  echo "${fpgm} :: ${rpgm} --> ${opgm}"
  pnmxarith -subtract -scale 3.0 -offset 0.5 ${fpgm} ${rpgm} > ${opgm}
}

rfile="reference.png"
frames=( `ls frame_[0-9][0-9][0-9][0-9][0-9].png | sort ` )
nf=${#frames[@]}     # Num frames.
cf=$(( ${nf} / 2 ))  # INdex of middle frame.
mfile="${frames[${cf}]}"

prepare_file ${rfile} .r.pgm 
prepare_file ${mfile} .m.pgm 
for ffile in "${frames[@]}" ; do
  orfile="${ffile/.png/_dr.pgm}"
  omfile="${ffile/.png/_dm.pgm}"
  echo "${ffile} --> ${orfile}, ${omfile}"
  prepare_file ${ffile} .f.pgm
  compare_files .f.pgm .r.pgm ${orfile}
  compare_files .f.pgm .m.pgm ${omfile}
  convert +append ${orfile} ${omfile} .drm.pgm
  display -title "${ffile}" -filter Box -resize '400%' .drm.pgm
done
