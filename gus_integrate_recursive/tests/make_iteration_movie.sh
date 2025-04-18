#! /bin/bash
# Last edited on 2025-04-12 19:17:31 by stolfi

outDir="$1"; shift
iterDir="$1"; shift
reportStep="$1"; shift
tag="$1"; shift
vmin="$1"; shift
vmax="$1"; shift

tmp="/tmp/$$"

echo "outDir = ${outDir}  iterDir = ${iterDir}" 1>&2 

files=( `( cd ${iterDir} && ls it-[0-9][0-9]-[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]-${tag}.fni ) | sort -t- -k2nr -k3n` )

echo "converting it-{level}-{iter}-${tag}.fni files to PNG ..." 1>&2
framePrefix="${tmp}-it"
nframes=0
for f in ${files[@]}; do
  # echo ${f} 1>&2
  levit="${f/-${tag}.fni/}"; levit="${levit/it-/}";
  level="${levit/-*/}"; xiter="${levit/*-/}";
  # echo "level = '${level}'  xiter = '${xiter}'" 1>&2
  niter=$(( 10#${xiter} + 0 ))
  rem=$(( ${niter} % ${reportStep} ))
  # echo "niter = '${niter}'  rem = '${rem}'" 1>&2
  if [[ ${rem} -eq 0 ]]; then
    fniFile="${iterDir}/${f}"
    xframe=`printf "%07d" "${nframes}"`
    plotFile="${framePrefix}-${xframe}.png"
    echo "  converting $f to ${plotFile} ..." 1>&2
    fni_plot.sh -channel 0 -range ${vmin} ${vmax} -title "level ${level} iter ${niter}" < ${fniFile} > ${plotFile}
    nframes=$(( ${nframes} + 1 ))
  fi
done

movieFile="${outDir}/${tag}-iters.mp4"
echo "creating movie ${movieFile} ..." 1>&2
ffmpeg \
  -framerate 10 \
  -start_number 0 \
  -i "${framePrefix}-%07d.png" \
  -y \
  -r 10 \
  -vcodec libx264 \
  -crf 25 \
  -pix_fmt yuv420p \
  ${movieFile}

mplayer ${movieFile}
