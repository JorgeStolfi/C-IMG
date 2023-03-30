#! /bin/bash
# Last edited on 2023-03-29 18:35:31 by stolfi

main="$1"; shift
name="$1"; shift

# POV-Ray paths:
POVRAY=/usr/bin/povray
POVINC=/usr/share/povray-3.7/include/
POVTTF=${STOLFIHOME}/tt-fonts
POVFIL=../povray_output_filter.gawk

function clean() {
  name="$1"; shift
  rm -f ${name}.png ${name}.log
}
        
function image() {
  main="$1"; shift
  name="$1"; shift
  width="$1"; shift
  height="$1"; shift
  
  echo "#declare scene_inc = \"${name}.inc\"" > params.inc 
  
  ln -s ../camera.inc
  ln -s ../isolines.ppm

  nrays=2
  clock=0.5000        
 
  ImageFile=${name}.png 

  rm -f ${name}.png
  nice ${POVRAY} \
      +FN +Q9 \
      +W${width} +H${height} \
      +AM1 +A0.0 +R${nrays} \
      +D \
      +K${clock} \
      +L${POVINC} \
      +L${POVTTF} \
      +I${main} \
      +O${name}.png \
      2>&1 \
    | ${POVFIL}
}

function showimage() {
  name="$1"; shift

  if [[ -s ${name}.png ]]; then
    display ${name}.png 
  fi
}

image ${main} ${name} 750 500
showimage ${name}
