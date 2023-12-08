#! /bin/bash 
# Last edited on 2023-08-27 21:53:43 by stolfi

# Reads a PGM or PBM image from {stdin}. Simulates the effect of the
# mapping from pixel coords to user coords or vice-versa, including axis
# flips and unit change but not origin shift. Writes the output image to
# {stdout}.

EXT="$1"; shift   # Image extension "ppm" or "pgm".
DIR="$1"; shift   # Direction: 0 input map, 1 output map
XAXIS="$1"; shift # X axis direction, "left" or "right".
YAXIS="$1"; shift # Y axis direction, "up" or "down".
UNIT="$1"; shift  # Pixels per user unit.

TYPE=`echo "${EXT}" | tr a-z A-Z` # Extension in uppercase.
ops=( )

# Simulate the "-xAxis" option:
if [[ "/${XAXIS}" == "/left" ]]; then
  ops+=( "|" "pnmflip" "-lr" )
elif [[ "/${XAXIS}" == "/right" ]]; then
  echo "no X flip" 1>& 2
else
  echo "** invalid XAXIS = '${XAXIS}'" 1>&2; exit 1
fi

# Simulate the "-yAxis" option:
if [[ "/${YAXIS}" == "/up" ]]; then
  ops+=( "|" "pnmflip" "-tb" )
elif [[ "/${YAXIS}" == "/down" ]]; then
  echo "no Y flip" 1>& 2
else
  echo "** invalid YAXIS = '${YAXIS}" 1>&2; exit 1
fi

if [[ ${DIR} -eq 0 ]]; then
  
  # Simulate the "-iUnit" opton:
  rpct=`echo "100.0/${UNIT}" |bc -lq`
  rpct=`echo "${rpct}" | sed -e 's:[.0][0]*$::g'`
  if [[ "/${rpct}" != "/1" ]]; then
    ops+=( "|" "convert" "${TYPE}:-" "-resize" "'${rpct}%'" "${TYPE}:-" )
  fi
elif [[ ${DIR} -eq 1 ]]; then  
  # Simulate the "-oUnit" opton:
  rpct=`echo "100.0*${UNIT}" | bc -lq`
  rpct=`echo "${rpct}" | sed -e 's:[.0][0]*$::g'`
  if [[ "/${rpct}" != "/1" ]]; then
    ops+=( "|" "convert" "${TYPE}:-" "-resize" "'${rpct}%'" "${TYPE}:-" )
  fi
else
  echo "** invalid DIR = '${DIR}" 1>&2; exit 1
fi

echo "commands = " "${ops[@]}" 1>&2

bash -c "cat ${ops[*]}"

