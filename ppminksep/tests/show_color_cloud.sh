#! /bin/bash
# Last edited on 2023-02-08 11:35:12 by stolfi

usage="$0 [-inGamma N.NN N.NN N.NN] [-logScale] IN_FILE OUT_PREFIX"

gawkopts=( )
outdir="out"
while [[ ($# -gt 0) && ("/$1" =~ /-*) ]]; do
  if [[ ( $# -ge 4 ) && ( "/$1" == "/-inGamma" ) ]]; then
    gawkopts=( ${gawkopts} -v "gmr=$2" -v "gmg=$3" -v "gmb=$4" )
    shift; shift; shift; shift
  elif [[ ( $# -ge 1 ) && ( "x$1" == "x-logScale" ) ]]; then
    gawkopts=( ${gawkopts} -v "logScale=1" )
    shift
  else
    echo "** unrecognized option $1" 1>&2; exit 1
  fi
done

if [[ $# -ne 2 ]]; then
  echo "** bad arguments line \"$@\"" 1>&2 
  echo "usage: ${usage}" 1>2; exit 1
fi

infile="$1"; shift;
outPrefix="$1"; shift;

ppmraw="${infile}"
parms="${outPrefix}.parms"

hist="${outPrefix}.hst"

gvpxs="${outPrefix}.gvpxs"
gvtet="${outPrefix}.gvtet"
gvtop="${outPrefix}.gvtop"

# Compute color histogram:
cat ${ppmraw} \
   | pnmscale 0.25 \
   | pnmdepth 127 \
   | pnmdepth 255 \
   | ppmhist \
   | gawk '/^ *[0-9]/{print;}' \
   | sort -b +3 -4nr \
   > ${hist}

# Create geomview file with pixel cloud:
cat ${hist} \
  | ./ppmhist_to_geomview.gawk \
      -f color_conv_funcs.gawk \
      ${gawkopts} \
  > ${gvpxs}

# Create geomview object showing color tetrahedron:
cat ${parms} \
  | ./parms_to_geomview.gawk \
      -f color_conv_funcs.gawk \
      ${gawkopts} \
  > ${gvtet}

# Create top-level geomview file:
/bin/rm -f ${gvtop}
echo "{ appearance {" >> ${gvtop}
echo "    lighting { ambient 1 1 1 }" >> ${gvtop}
echo "  }" >> ${gvtop}
echo "  LIST" >> ${gvtop}
echo "  { geom { < ${gvpxs} } } " >> ${gvtop}
echo "  { geom { < ${gvtet} } } " >> ${gvtop}
echo "}" >> ${gvtop}

# Run geomview:
geomview -b 0.9 0.9 0.9 -c '(bbox-draw allgeoms no)' ${gvtop}

# Cleanup:
rm ${gvpxs} ${gvtet} ${gvtop}
