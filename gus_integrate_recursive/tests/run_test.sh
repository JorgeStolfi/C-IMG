#! /bin/bash
# Last edited on 2025-04-08 09:11:29 by stolfi

testName="$1"; shift;        # E.g. "cbabel".
size="$1"; shift;            # E.g. "0256x0192"
hintsWeight="$1"; shift;     # In {[0 _ 1]}. 
hasRefZ="$1"; shift;         # 0 or 1
initial_option="$1"; shift;  # "00", "MF", or "SY"
initial_noise="$1"; shift;   # Fraction in 0 _ 1.
maxLevel="$1"; shift         # 0 or more.

progDir=".."
prog="gus_integrate_recursive"

inDir="in/${testName}/${size}"

outDir="out/${testName}/${size}"
mkdir -p ${outDir}

# Runs ${prog} on one set of images.
#
# Inputs, in directory ${inDir}:
#
#   ${size}-PM-G.fni  Gradient (slope) map from photometric stereo.
#   ${size}-PM-N.fni  Normal map from photometric stereo.
#   ${size}-MF-Z.fni  Height map from multifocus stereo.
#   ${size}-RF-Z.fni  Reference (ideal) height map for evaluation.
#
# The normal map is used only if the gradient map is not present.
# The multifocus and reference height maps are optional.
#
# Outputs, in ${outDir}/${testName}:
#
#   ${size}-E.fni
#   ???

tmp="/tmp/$$"
mkdir -p ${outDir}

in_slopes_fni="${inDir}/PM-G.fni"
in_normals_fni="${inDir}/PM-N.fni"
if [[ -s ${in_slopes_fni} ]]; then
  input_options=( -slopes ${in_slopes_fni} )
elif [[ -s ${in_normals_fni} ]]; then
  input_options=( -normals ${in_normals_fni} )
else
  echo "** ${in_slopes_fni} and ${in_normals_fni} are both missing" 1>&2 ; exit 1
fi

if [[ ${hasRefZ} -ne 0 ]]; then
  in_refz_fni="${inDir}/RF-Z.fni"
  out_errz_fni="${outPrefix}-00-end-E.fni"
  reference_options=( "-reference" ${in_refz_fni} )
else
  in_refz_fni="NONE/NONE"
  out_errz_fni="NONE/NONE"
  reference_options=( )
fi

zeroPat=':[0.]*:'
if [[ ! ( ":${hintsWeight}:" =~ ${zeroPat} ) ]]; then
  in_hints_fni="${inDir}/MF-Z.fni"
  hints_options=( "-hints" ${in_hints_fni} ${hintsWeight} )
else
  hints_options=()
fi

if [[ "/${initial_option}" == "/00" ]]; then
  initial_option=( "zero" )
elif [[ "/${initial_option}" == "/MF" ]]; then
  initial_option=( "hints" )
elif [[ "/${initial_option}" == "/RF" ]]; then
  initial_option=( "reference" )
else
  echo "** invalid initial_option = '${initial_option}'" 1>&2 ; exit 1
fi

if [[ ${size} == "0512x0422" ]]; then
  clear_options=( -clear 0 511 0 40 )
else
  clear_options=( )
fi

outPrefix="${outDir}/out"
iterPrefix="${outDir}-iters/it"

rm -f {${outPrefix},${iterPrefix}}*.{fni,sys,pgm,txt,png}
set -x
${progDir}/${prog} \
  ${input_options[@]} \
  ${hints_options[@]}  \
  ${reference_options[@]} \
  ${clear_options[@]} \
  -maxLevel ${maxLevel} \
  -initial ${initial_option} ${initial_noise} \
  -outPrefix ${outPrefix} \
  -reportStep 10 \
  -verbose
set +x

function showfni() {
  tgfi="$1"; shift

  tag="${tgfi/:*}"
  fni_file="${tgfi/*:}"
  fni_name="${fni_file/*\//}"
  
  if [[ -s ${fni_file} ]]; then
    echo "showing ${fni_file} ..." 1>&2
    ny="`cat ${fni_file} | egrep -e 'NY *='`"
    ny="${ny/*=/}"; ny="${ny/* /}"; 
    ny=$(( 10#${ny} + 0 ))

    if [[ "/${tag/* /}" == "/G"  ]]; then
      scale=1;
    else 
      scale=$(( ${ny} / 2 ));
    fi
    echo "tag = '${tag}'  ny = ${ny}  scale=${scale}" 1>&2
    pgm_file="${tmp}.pgm"
    png_file="${fni_file/.fni/.png}"
    fni_view -title "${tag} : ${fni_name}" -hist y -scale ${scale} ${fni_file}
    # fni_to_pnm -channel 0 -yAxis up < ${fni_file} > ${pgm_file}
    # convert ${pgm_file} ${png_file}
    # display -title "${tag} : ${fni_name}" -filter box -resize 'x800<' ${png_file}
  else
    echo "** ${fni_file} not present" 1>&2
  fi
}

for tgfi in 'input PM-G':${in_slopes_fni} 'input PM-N':${in_normals_fni} 'input MF-Z':${in_hints_fni} 'input RF-Z':${in_refz_fni} ; do
  showfni "${tgfi}"
done

for level in 20 19 18 17 16 15 14 13 12 11 10 09 08 07 06 05 04 03 02 01 00 ; do
  levelPrefix="${outPrefix}-${level}"
  for state in "beg" "end"; do 
    for tag in 'Z' 'G' 'E' ; do
      if [[ ( "/${state}" == "/beg" ) || ( "/${tag}" != "/G" ) ]]; then
        fni_file="${levelPrefix}-${state}-${tag}.fni"
        stg="${state} ${tag}"
        showfni "${stg}:${fni_file}"
      fi
    done
  done
done

make_iteration_movie.sh ${outPrefix} 10 Z -${zscale} +${zscale}
make_iteration_movie.sh ${outPrefix} 10 E -2.0 +2.0
echo "OK"

